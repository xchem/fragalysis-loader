FROM xchem/fragalysis-backend
ENV PYTHONUNBUFFERED 1
ADD . /code/
WORKDIR /code
