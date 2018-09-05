#!groovy​

// The 'xchem' Fragalysis Loader Jenkinsfile.

pipeline {

  agent { label 'buildah-slave' }

  environment {
    // Local registry details (for FROM image)
    REGISTRY_USER = 'jenkins'
    REGISTRY = 'docker-registry.default:5000'
    STREAM_IMAGE = "${REGISTRY}/fragalysis-cicd/fragalysis-loader:latest"

    // Slack channel for all notifications
    SLACK_BUILD_CHANNEL = 'dls-builds'
    // Slack channel to be used for errors/failures
    SLACK_ALERT_CHANNEL = 'dls-alerts'
  }

  stages {

    stage('Build Image') {
      steps {
        slackSend channel: "#${SLACK_BUILD_CHANNEL}",
                  message: "${JOB_NAME} build ${BUILD_NUMBER} - starting..."
        script {
          TOKEN = sh(script: 'oc whoami -t', returnStdout: true).trim()
        }
        sh "buildah bud --tls-verify=false --creds=${REGISTRY_USER}:${TOKEN} --format docker -f Dockerfile-cicd -t ${STREAM_IMAGE} ."
      }
    }

    stage('Push Image') {
      steps {
        script {
          TOKEN = sh(script: 'oc whoami -t', returnStdout: true).trim()
        }
        sh "podman login --tls-verify=false --username ${REGISTRY_USER} --password ${TOKEN} ${REGISTRY}"
        sh "buildah push --tls-verify=false ${STREAM_IMAGE} docker://${STREAM_IMAGE}"
        sh "podman logout ${REGISTRY}"
      }
    }

  }

  // Post-job actions.
  // See https://jenkins.io/doc/book/pipeline/syntax/#post
  post {

    success {
      slackSend channel: "#${SLACK_BUILD_CHANNEL}",
                color: 'good',
                message: "${JOB_NAME} build ${BUILD_NUMBER} - complete"
    }

    failure {
      slackSend channel: "#${SLACK_BUILD_CHANNEL}",
                color: 'danger',
                message: "${JOB_NAME} build ${env.BUILD_NUMBER} - failed (${BUILD_URL})"
      slackSend channel: "#${SLACK_ALERT_CHANNEL}",
                color: 'danger',
                message: "${JOB_NAME} build ${BUILD_NUMBER} - failed (${BUILD_URL})"
    }

    fixed {
      slackSend channel: "#${env.SLACK_ALERT_CHANNEL}",
                color: 'good',
                message: "${JOB_NAME} build - fixed"
    }

  }

}
