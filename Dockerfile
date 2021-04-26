ARG BE_NAMESPACE="xchem"
ARG BE_IMAGE_TAG="latest"
FROM ${BE_NAMESPACE}/fragalysis-backend:${BE_IMAGE_TAG}
ENV PYTHONUNBUFFERED 1

# Build provenance.
# Add labels for the build arguments...
LABEL BE_NAMESPACE=${BE_NAMESPACE} \
      BE_IMAGE_TAG=${BE_IMAGE_TAG}

ADD . /code/
WORKDIR /code
