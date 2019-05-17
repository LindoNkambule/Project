#First pull the docker2singularity image using the command below:
#docker run singularityware/docker2singularity

docker run -v /var/run/docker.sock:/var/run/docker.sock \
-v /path/to/directory:/output \
--privileged -t --rm \
singularityware/docker2singularity \
--name honours_project ubuntu:18.04

#-v specify the path in which you want to save the singularity image to
#--name specify the name you want to give to the image
