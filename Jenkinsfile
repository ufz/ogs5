#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.15') _

pipeline {
    agent none
    stages {
        stage('Parallel Stage') {
            parallel {
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
                                tools: [msBuild(name: 'MSVC', pattern: 'build/build.log')],
                                unstableTotalAll: 1
                        }
                        success {
                            archiveArtifacts 'build/*.zip'
                        }
                    }
                }
                stage('Linux') {
                    agent { label 'envinf1' }
                    steps {
                        dir('benchmarks') { git 'https://github.com/ufz/ogs5-benchmarks.git' }
                        dir('benchmarks_ref') { git 'https://github.com/ufz/ogs5-benchmarks_ref.git' }
                        script {
                            configure {
                                cmakeOptions =
                                    "-DOGS_CONFIG=FEM " +
                                    "-DOGS_USE_CVODE=ON " +
                                    "-DOGS_NO_EXTERNAL_LIBS=ON " +
                                    "-DNUMDIFF_TOOL_PATH=/usr/local/numdiff/5.8.1-1/bin/numdiff" +
                                    "-DOGS_CPU_ARCHITECTURE=generic"
                                dir="build_FEM"
                                env="envinf1/cli.sh"
                            }
                            build {
                                dir="build_FEM"
                                env="envinf1/cli.sh"
                                log="build.log"
                            }
                            build {
                                dir="build_FEM"
                                env="envinf1/cli.sh"
                                target="benchmarks-short-normal-long"
                            }
                        }
                    }
                    post {
                        always {
                            recordIssues enabledForFailure: true,
                                tools: [gcc4(name: 'GCC', pattern: 'build_FEM/build.log')],
                                unstableTotalAll: 3
                            xunit([CTest(pattern: 'build_FEM/Testing/**/*.xml')])
                            archiveArtifacts allowEmptyArchive: true,
                                             artifacts: 'build_FEM/benchmarks/**/*.numdiff'
                        }
                        success {
                            archiveArtifacts 'build_FEM/*.tar.gz'
                        }
                    }
                }
            }
        }
    }
}
