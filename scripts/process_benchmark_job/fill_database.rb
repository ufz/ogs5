require 'rubygems'
require 'sequel'
$DB = Sequel.connect('sqlite://' + ARGV[0])

require 'commit_info.rb'
require 'bench_info.rb'


$DB.create_table! :commit_infos do
  primary_key     :revision
  Time            :date
  String          :branch
  foreign_key     :author_id, :table => :authors
end

$DB.create_table! :benchmark_runs do
  primary_key :id
  String      :name
  Float       :time
  Boolean     :crashed
  Boolean     :passed
  String      :config

  foreign_key :author_id
  index       :author_id
  foreign_key :commit_info_id
  index       :commit_info_id
end

commit = CommitInfoLoader.new('tests/svnInfoOld.txt')
info = BenchmarkRunsLoader.new('tests/benchOutOld.txt')
commit = CommitInfoLoader.new('tests/svnInfo.txt')
info = BenchmarkRunsLoader.new('tests/benchOut.txt')

$DB[:benchmark_runs].each {|row| p row}
$DB[:commit_infos].each {|row| p row}