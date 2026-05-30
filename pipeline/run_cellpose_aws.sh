#!/usr/bin/env bash
set -euo pipefail
REGION="us-east-2"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

EC2_TASK_ROLE="cellpose-task-role"
EC2_INSTANCE_PROFILE="cellpose-profile"

S3_BUCKET_MODEL="stan-img-processing"
S3_BUCKET_DATA="stan-data-misc"
EC2_TRUST_POLICY=$(cat <<'POLICY'
{
    "Version": "2012-10-17",
    "Statement": [{
      "Effect": "Allow",
      "Principal": {"Service": "ec2.amazonaws.com"},
      "Action": "sts:AssumeRole"
    }]
}
POLICY
)

S3_POLICY=$(cat <<'POLICY'
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "s3:GetObject",
                "s3:PutObject",
                "s3:ListBucket",
                "s3:DeleteObject"
            ],
            "Resource": [
                "arn:aws:s3:::stan-img-processing",
                "arn:aws:s3:::stan-img-processing/*",
                "arn:aws:s3:::stan-data-misc",
                "arn:aws:s3:::stan-data-misc/*"
            ]
        }
    ]
}
POLICY
)

CELLPOSE_ENV=$(cat <<'ENV'
--extra-index-url https://download.pytorch.org/whl/cu126
cellpose==3.1.0
numpy==2.0.2
scipy==1.15.2
opencv-python-headless==4.11.0.86
tifffile==2025.2.18
imagecodecs==2024.12.30
fastremap==1.15.1
numba==0.61.0
llvmlite==0.44.0
torch==2.6.0+cu126
torchvision==0.21.0+cu126
scikit-learn==1.5.2
natsort==8.4.0
tqdm==4.67.1
roifile==2025.2.20
pillow==11.0.0
matplotlib==3.10.0
ENV
)

create_role() {
    local role_name="$1"
    local trust_policy="$2"
    if ! aws iam get-role --role-name "${role_name}" &>/dev/null; then
        echo "Create role ${role_name}"
        aws iam create-role \
            --role-name "${role_name}" \
            --assume-role-policy-document "${trust_policy}"
    else
        echo "Role ${role_name} exists. Updating role policy"
        aws iam update-assume-role-policy \
            --role-name "${role_name}" \
            --assume-role-policy-document "${trust_policy}"
    fi
    echo "Done"
}

create_role "${EC2_TASK_ROLE}" "${EC2_TRUST_POLICY}"
aws iam put-role-policy \
    --role-name "${EC2_TASK_ROLE}" \
    --policy-name "s3-access" \
    --policy-document "${S3_POLICY}"

aws iam create-instance-profile \
    --instance-profile-name $EC2_INSTANCE_PROFILE

aws iam add-role-to-instance-profile \
    --instance-profile-name $EC2_INSTANCE_PROFILE \
    --role-name $EC2_TASK_ROLE

CELLPOSE_MODEL="olympus_scope_dapi_20x_epoch_2800"
DIAMETER=26
DO_3D=true
ANISOTROPY=6
CELLPROB_THRESHOLD=-1
FLOW_THRESHOLD=0.5
PLATE="plate1_1"

AMI_ID=$(aws ec2 describe-images --region $REGION \
    --filters "Name=name,Values=*Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu*" \
              "Name=state,Values=available" \
    --query 'Images | sort_by(@, &CreationDate) | [-1].ImageId' \
    --output text)
echo "Using AMI: $AMI_ID"

KEY_NAME="cellpose-key"
if ! aws ec2 describe-key-pairs --region $REGION --key-names $KEY_NAME &>/dev/null; then
    aws ec2 create-key-pair --region $REGION --key-name $KEY_NAME \
        --query 'KeyMaterial' --output text > ~/.ssh/${KEY_NAME}.pem
    chmod 600 ~/.ssh/${KEY_NAME}.pem
    echo "Created key pair: ~/.ssh/${KEY_NAME}.pem"
fi

SG_NAME="cellpose-sg"
SG_ID=$(aws ec2 describe-security-groups --region $REGION \
    --filters "Name=group-name,Values=$SG_NAME" \
    --query 'SecurityGroups[0].GroupId' --output text 2>/dev/null)
if [ "$SG_ID" = "None" ] || [ -z "$SG_ID" ]; then
    SG_ID=$(aws ec2 create-security-group --region $REGION \
        --group-name $SG_NAME \
        --description "Cellpose spot instance" \
        --query 'GroupId' --output text)
    aws ec2 authorize-security-group-ingress --region $REGION \
        --group-id $SG_ID --protocol tcp --port 22 --cidr 0.0.0.0/0
    echo "Created security group: $SG_ID"
fi

echo "$CELLPOSE_ENV" > /tmp/cellpose_requirements.txt
aws s3 cp /tmp/cellpose_requirements.txt s3://$S3_BUCKET_MODEL/cellpose-models/cellpose_requirements.txt

USER_DATA=$(cat <<'USERDATA'
#!/bin/bash
set -euxo pipefail
exec > /var/log/cellpose.log 2>&1

WORK=/opt/cellpose
mkdir -p $WORK && cd $WORK

# wait for GPU driver
for i in $(seq 1 30); do nvidia-smi && break || sleep 10; done

# create venv and install
python3 -m venv venv
source venv/bin/activate
aws s3 cp s3://stan-img-processing/cellpose-models/cellpose_requirements.txt requirements.txt
pip install -r requirements.txt

# download model
mkdir -p /root/.cellpose/models
aws s3 cp s3://stan-img-processing/cellpose-models/olympus_scope_dapi_20x_epoch_2800 \
    /root/.cellpose/models/olympus_scope_dapi_20x_epoch_2800

# download data
mkdir -p $WORK/data
aws s3 sync s3://stan-data-misc/cellpose-2026-02-24/plate1_1/ $WORK/data/

# run cellpose
cellpose --dir $WORK/data \
    --pretrained_model olympus_scope_dapi_20x_epoch_2800 \
    --diameter 26 \
    --do_3D --anisotropy 6 \
    --cellprob_threshold -1 \
    --flow_threshold 0.5 \
    --use_gpu \
    --save_tif \
    --verbose

# upload results
aws s3 sync $WORK/data/ s3://stan-data-misc/cellpose-2026-02-24/plate1_1/ \
    --exclude "*.tif" \
    --include "*_cp_masks*" \
    --include "*_seg.npy"

# shutdown
shutdown -h now
USERDATA
)

sleep 10

INSTANCE_ID=$(aws ec2 run-instances --region $REGION \
    --image-id $AMI_ID \
    --instance-type g4dn.4xlarge \
    --instance-market-options '{"MarketType":"spot"}' \
    --iam-instance-profile "Name=$EC2_INSTANCE_PROFILE" \
    --key-name $KEY_NAME \
    --security-group-ids $SG_ID \
    --user-data "$USER_DATA" \
    --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":100,"VolumeType":"gp3"}}]' \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=cellpose-${PLATE}}]" \
    --query 'Instances[0].InstanceId' \
    --output text)

echo "Launched spot instance: $INSTANCE_ID"
echo "Monitor: aws ec2 describe-instances --region $REGION --instance-ids $INSTANCE_ID --query 'Reservations[0].Instances[0].State.Name'"
echo "SSH: ssh -i ~/.ssh/${KEY_NAME}.pem ubuntu@\$(aws ec2 describe-instances --region $REGION --instance-ids $INSTANCE_ID --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)"
echo "Logs on instance: sudo cat /var/log/cellpose.log"
