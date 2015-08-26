def flow
node('envinf1')
{
	dir('ogs') {
		checkout scm
	}

	checkout 'https://github.com/ufz/ogs5-benchmarks.git', 'master', 'benchmarks'
	checkout 'https://github.com/ufz/ogs5-benchmarks_ref.git', 'master', 'benchmarks_ref'

	flow = load 'ogs/scripts/jenkins/linux.groovy'
	flow.runLinux()

} // end node envinf1

// Check out branch from url into directory
def checkout(url, branch, directory) {
	stage 'Checkout'
	checkout([$class: 'GitSCM',
		branches: [[name: "*/${branch}"]],
		doGenerateSubmoduleConfigurations: false,
		extensions:
			[[$class: 'RelativeTargetDirectory', relativeTargetDir: "${directory}"]],
		submoduleCfg: [],
		userRemoteConfigs:
			[[credentialsId: '6c1dad0d-9b3c-44c2-a6bb-669562045187', url: "${url}"]]])
}
