#First pull the docker2singularity image using the command below:
#docker run singularityware/docker2singularity

docker run -v /var/run/docker.sock:/var/run/docker.sock \
-v /Users/lindokuhle/Desktop/Project:/output \
--privileged -t --rm \
singularityware/docker2singularity \
--name honours_project ubuntu:18.04
