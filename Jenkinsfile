#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.15') _

pipeline {
    agent none
    stages {
        stage('Win') {
            agent { label 'win' }
            environment {
                MSVC_NUMBER = '15'
                MSVC_VERSION = '2017'
            }
            steps {
                script {
                    configure { }
                    build { log="build.log" }
                }
            }
            post {
                always {
                    recordIssues enabledForFailure: true, filters: [
                        excludeFile('.*\\.conan.*'), excludeFile('.*ThirdParty.*'),
                        excludeFile('.*thread.hpp')],
                        tools: [msBuild(name: 'MSVC', pattern: 'build/build.log')]
                }
                success {
                    archiveArtifacts 'build/*.zip'
                }
            }

        }
    }
}
