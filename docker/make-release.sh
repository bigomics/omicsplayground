## run from root

## get version
version=`head -n1 VERSION`
image=bigomics/omicsplayground
echo VERSION=$version

## update datasets info
echo '(cd scripts && Rscript update-datasets-info.R)'

##---------------- github ----------------------

echo git pull
echo git tag -a \"$version\" -m \"version $version\"
echo git push
echo git push --tags

##---------------- docker ----------------------
echo "nohup docker build -f docker/Dockerfile.base --no-cache -t $image:base . > docker.out &"
echo "nohup docker build -f docker/Dockerfile.testing --no-cache -t $image:testing . > docker.out &"
echo "nohup docker build -f docker/Dockerfile.testing2 --no-cache -t $image:testing . > docker.out &"
echo "nohup docker build -f omicsplayground/docker/Dockerfile.dev --no-cache -t $image:dev . > docker.out &"

## +run in background, remove contained after use, give nice name
echo docker run --rm -d -p 4000:3838 --name=play1 $image

## enter bash in the container
echo docker exec -it play1 /bin/bash
echo docker run -it --entrypoint /bin/bash bigomics/omicsplayground:testing -s

## stop the container
echo docker stop play1

## To save/load your Docker Image as a tar-archive
echo docker save -o ~/omicsplayground_$version.tar $image
echo docker load -i ~/omicsplayground_$version.tar

## give version tag
echo docker tag $image:latest $image:$version

## publish to Docker Hub
echo docker login
echo docker push bigomics/omicsplayground:$version
echo docker push bigomics/omicsplayground:latest
