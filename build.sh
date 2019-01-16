set -ex
# SET THE FOLLOWING VARIABLES
# docker hub username
USERNAME=lindonkambule
# image name
IMAGE=project
docker build -t $USERNAME/$IMAGE:latest .
