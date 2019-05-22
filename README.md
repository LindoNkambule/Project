# Multiple OpenJDK versions
The Docker image has multiple java versions: JDK 8 and JDK 11. This is because JDK 11 is the default version for Ubuntu 18.04 and GATK runs on java version 1.8.x (it is not compatible with other java versions).

To check the versions, inside the image run:
```
update-java-alternatives --list
```

If you want to use only GATK, you have to set JDK 8 as the default version. Run the following command:
```
update-alternatives --config java
#Select option 2, JDK 8.
```



If you want to run a variant calling pipeline, from data preparation, alignment, to variant calling, please set a java variable (for JDK 11) as shown in line 12 of the 'whole_pipeline.sh' script. This is because PICARD tools gives an error when you run it using JDK 8. I am still working on finding what's causing the issue.


# Running variant calling with GATK using nextflow and Docker
There is an image on Docker Hub with GATK already installed. To run the variant calling pipeline with GATK using nextflow and Docker, using the pipeline.nf script, use the command:
```
nextflow run pipeline.nf -with-docker lindonkambule/project_gatk
```

# Running variant calling with GATK using nextflow and Singularity
With the command below, nextflow will run the variant calling pipeline using Singularity but with the specified Docker image. Yes, Docker image. The 'docker://' option will download the specified Docker image from Docker Hub. This image will then be converted to a Singularity image for the -with-singularity option to work.
```
nextflow run pipeline.nf -with-singularity docker://lindonkambule/project_gatk
```
