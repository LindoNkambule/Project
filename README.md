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


# Running variant calling using GATK and nextflow
There is a an image on dockerhub with GATK already installed. To run the variant calling pipeline using GATK and nextflow, using the pipeline.nf script, use the command:
```
nextflow run pipeline.nf -with-docker lindonkambule/project_gatk
```
