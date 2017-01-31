defaultCMakeOptions = '-DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=$CMAKE_TOOLCHAIN_FILE'

node('docker') {
    stage('Checkout') {
        dir('ogs') { checkoutWithTags('https://github.com/ufz/ogs5.git') }
    }

    docker.image('ogs6/mingw-base').inside() {
        stage('Configure') {
            configure 'build', ''
        }

        stage('Build') {
            build 'build', ''
            if (env.BRANCH_NAME == 'master')
                build 'build', 'package'
        }
    }

    stage('Post') {
        archive 'build*/*.zip'
    }
}

def configure(buildDir, cmakeOptions) {
    sh "rm -rf ${buildDir} && mkdir ${buildDir}"
    sh "cd ${buildDir} && cmake ../ogs ${defaultCMakeOptions} ${cmakeOptions}"
}

def build(buildDir, target) {
    sh "cd ${buildDir} && make -j \$(nproc) ${target}"
}
