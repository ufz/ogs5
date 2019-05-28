#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.16') _

pipeline {
    agent none
    options {
      ansiColor('xterm')
      timestamps()
      buildDiscarder(logRotator(numToKeepStr: '30', artifactNumToKeepStr: '10'))
      timeout(time: 3, unit: 'HOURS')
    }
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
                stage('Linux-FEM') {
                    agent { label 'frontend2' }
                    steps {
                        script {
                            ogs5BuildLinux {
                                config="FEM"
                                cmakeOptions = " -DOGS_USE_CVODE=ON " +
                                               " -DOGS_NO_EXTERNAL_LIBS=ON "
                            }
                        }
                    }
                    post {
                        always { script { ogs5PostAlways { config="FEM" } } }
                        success { script { ogs5PostSuccess { config="FEM" } } }
                    }
                }
                stage('Linux-SP') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux { config="SP" } }
                    }
                    post {
                        always { script { ogs5PostAlways { config="SP" } } }
                        success { script { ogs5PostSuccess { config="SP" } } }
                    }
                }
                stage('Linux-GEMS') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux { config="GEMS" } }
                    }
                    post {
                        always { script { ogs5PostAlways {
                            config="GEMS"
                            warnings=36 } } }
                        success { script { ogs5PostSuccess { config="GEMS" } } }
                    }
                }
                stage('Linux-PQC') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux { config="PQC" } }
                    }
                    post {
                        always { script { ogs5PostAlways { config="PQC" } } }
                        success { script { ogs5PostSuccess { config="PQC" } } }
                    }
                }
                stage('Linux-IPQC') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux { config="IPQC" } }
                    }
                    post {
                        always { script { ogs5PostAlways { config="IPQC" } } }
                        success { script { ogs5PostSuccess { config="IPQC" } } }
                    }
                }
                stage('Linux-BRNS') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux { config="BRNS" } }
                    }
                    post {
                        always { script { ogs5PostAlways {
                            config="BRNS"
                            warnings=3 } } }
                        success { script { ogs5PostSuccess { config="BRNS" } } }
                    }
                }
                stage('Linux-MKL') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux {
                            config="MKL"
                            cmakeOptions="-DMKL_DIR=/opt/intel/mkl"
                        } }
                    }
                    post {
                        always { script { ogs5PostAlways { config="MKL" } } }
                        success { script { ogs5PostSuccess { config="MKL" } } }
                    }
                }
                stage('Linux-LIS') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux { config="LIS" } }
                    }
                    post {
                        always { script { ogs5PostAlways {
                            config="LIS"
                            warnings=0 } } }
                        success { script { ogs5PostSuccess { config="LIS" } } }
                    }
                }
                stage('Linux-MPI') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux { config="MPI" } }
                    }
                    post {
                        always { script { ogs5PostAlways { config="MPI" } } }
                        success { script { ogs5PostSuccess { config="MPI" } } }
                    }
                }
                stage('Linux-PETSC') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux { config="PETSC" } }
                    }
                    post {
                        always { script { ogs5PostAlways {
                            config="PETSC"
                            warnings=4 } } }
                        success { script { ogs5PostSuccess { config="PETSC" } } }
                    }
                }
                stage('Linux-PETSC_GEMS') {
                    agent { label 'frontend2' }
                    steps {
                        script { ogs5BuildLinux { config="PETSC_GEMS" } }
                    }
                    post {
                        always { script { ogs5PostAlways {
                            config="PETSC_GEMS"
                            warnings=36 } } }
                        success { script { ogs5PostSuccess { config="PETSC_GEMS" } } }
                    }
                }
            }
        }
    }
}
