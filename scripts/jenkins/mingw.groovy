defaultCMakeOptions = '-DCMAKE_TOOLCHAIN_FILE=$CMAKE_TOOLCHAIN_FILE'

node('docker')
{
    step([$class: 'GitHubSetCommitStatusBuilder', statusMessage: [content: 'Started Jenkins mingw build']])

    stage 'Checkout'
    dir('ogs') { checkout scm }

    docker.image('ogs6/mingw-base').inside()
    {
        stage 'Configure'
        configure 'build', ''

        stage 'Build'
        build 'build', ''
        // if (env.BRANCH_NAME == 'master')
        build 'build', 'package'
    }

    stage 'Post'
    archive 'build*/*.zip'
    step([$class: 'S3BucketPublisher', dontWaitForConcurrentBuildCompletion: false, entries: [[bucket: 'opengeosys/ogs5-binaries/head', excludedFile: '', flatten: true, gzipFiles: false, managedArtifacts: false, noUploadOnFailure: true, selectedRegion: 'eu-central-1', sourceFile: 'build*/*.zip', storageClass: 'STANDARD', uploadFromSlave: false, useServerSideEncryption: false]], profileName: 'S3 UFZ', userMetadata: []])


    step([$class: 'GitHubCommitNotifier', resultOnFailure: 'FAILURE', statusMessage: [content: 'Finished Jenkins mingw build']])
}

def configure(buildDir, cmakeOptions) {
    sh "rm -rf ${buildDir} && mkdir ${buildDir}"
    sh "cd ${buildDir} && cmake ../ogs ${defaultCMakeOptions} ${cmakeOptions}"
}

def build(buildDir, target) {
    sh "cd ${buildDir} && make -j \$(nproc) ${target}"
}

properties ([[$class: 'BuildDiscarderProperty', strategy: [$class: 'LogRotator', artifactDaysToKeepStr: '', artifactNumToKeepStr: '5', daysToKeepStr: '', numToKeepStr: '25']]])
