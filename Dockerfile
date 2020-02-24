ARG from_image="xchem/fragalysis-backend"
ARG from_tag="latest"
FROM ${from_image}:${from_tag}
ENV PYTHONUNBUFFERED 1
ADD . /code/
WORKDIR /code
