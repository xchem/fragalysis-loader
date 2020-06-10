ARG BE_NAMESPACE="xchem"
ARG BE_IMAGE_TAG="latest"
FROM ${BE_NAMESPACE}/fragalysis-loader:${BE_IMAGE_TAG}

COPY docker-entrypoint.sh ./
COPY loader-requirements.txt ./
RUN pip install -r loader-requirements.txt && \
    chmod -R 755 *.sh && \
    chmod -R 755 *.py

CMD ["./docker-entrypoint.sh"]
