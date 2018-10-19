def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()

stage('Checkout') {
    dir('benchmarks') { git 'https://github.com/ufz/ogs5-benchmarks.git' }
    dir('benchmarks_ref') { git 'https://github.com/ufz/ogs5-benchmarks_ref.git' }
}

def configs = [
    [name: "FEM", cmakeOptions: " -DOGS_USE_CVODE=ON -DOGS_NO_EXTERNAL_LIBS=ON ",
        target: "package", artifacts: "*.tar.gz"],
    [name: "SP"],
    [name: "GEMS"],
    [name: "PQC"],
    [name: "IPQC"],
    [name: "BRNS"],
    [name: "MKL", cmakeOptions: " -DMKL_DIR=/opt/intel/mkl "],
    [name: "LIS"],
    [name: "MPI"],
    [name: "PETSC"],
    [name: "PETSC_GEMS"]
]

stage('Build Configs') {
    def buildTasks = [:]
    for (i = 0; i < configs.size(); i++) {
        def config = configs[i]
        def defaultCMakeOptions =
            ' -DOGS_CONFIG=' + config.name +
            ' -DCMAKE_BUILD_TYPE=Release' +
            ' -DNUMDIFF_TOOL_PATH=/usr/local/numdiff/5.8.1-1/bin/numdiff' +
            ' -DOGS_CPU_ARCHITECTURE=generic'
        def cmakeOptions = defaultCMakeOptions +
            (config.cmakeOptions ? config.cmakeOptions : '')

        buildTasks[config.name] = {
            configure.linux(
                cmakeOptions: cmakeOptions,
                dir: 'build_' + config.name,
                env: 'envinf1/cli.sh',
                script: this)
            build.linux(
                cmd: 'make -j 4',
                dir: 'build_' + config.name,
                env: 'envinf1/cli.sh',
                script: this,
                target: config.target ? config.target : '')

            if (config.artifacts && env.BRANCH_NAME == 'master')
                archive 'build_' + config.name + '/' + config.artifacts
        }
    }
    parallel buildTasks
}

stage('Benchmarking') {
    def benchmarkTasks = [:]
    for (i = 0; i < configs.size(); i++) {
        def config = configs[i]
        benchmarkTasks[config.name] = {
            catchError {
                build.linux(
                    cmd: 'make -j 4',
                    dir: 'build_' + config.name,
                    env: 'envinf1/cli.sh',
                    script: this,
                    target: 'benchmarks-short-normal-long')
            }
        }
    }
    parallel benchmarkTasks
    archive '**/*.numdiff'
}

stage('Post') {
    step([$class: 'LogParserPublisher', failBuildOnError: true, unstableOnWarning: true,
        projectRulePath: 'ogs/scripts/jenkins/log-parser.rules', useProjectRule: true])
    for (i = 0; i < configs.size(); i++) {
        sh("""shopt -s extglob
              rm -rf build_${configs[i].name}/!(benchmarks)
        """.stripIndent())
    }
}
