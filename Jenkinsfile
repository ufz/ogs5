node('master') {
    checkout scm
    parallel linux: { load 'scripts/jenkins/linux.groovy' },
    mingw: { load 'scripts/jenkins/mingw.groovy' }
}
