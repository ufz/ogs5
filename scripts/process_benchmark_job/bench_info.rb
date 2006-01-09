require 'rubygems'
require 'sequel'
require 'time'
#require 'database.rb'
require './ogs_author_mapping.rb'
require './commit_info.rb'

$DB.create_table? :benchmark_runs do
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

class BenchmarkRun < Sequel::Model( :benchmark_runs )
  many_to_one :author
  many_to_one :commit_info

  def benchmark_info
    return "#{self.name} by #{self.author.name}"
  end

  def inspect2
    if @crashed
      passed = 'crashed'
    elsif @passed
      passed = 'passed'
    else
      passed = 'failed'
    end
    @passed ? passed = 'passed' : passed = 'failed'
    if not @passed
      p self if not self.name
      p self if not self.author
      puts "Benchmark #{benchmark_info} #{passed} in #{time} s."
    end
  end
end


class BenchmarkRunsLoader

  attr_reader :bench_test_infos

  def initialize(filename, commit_info)
    @bench_test_infos = []
    puts "Loading benchmark job from #{filename}"
    File.open(filename, 'r') do |file|
      num_test_project_lines = 0
      config = nil
      while line = file.gets
        if line =~ /Test project/
          # Check config
          config = line.scan(/build_([\w]+$)/)[0].to_s

          # Even test runs are benchmarks, otherwise file compares
          num_test_project_lines += 1

        elsif line =~ /^\s*[0-9]+\/[0-9]+\sTest/
          # Check for passed or failed
          crashed = false
          crashed = false if line =~ /\s+Passed\s+/
          crashed = true if line =~ /\*+Failed\s+/

          # Check author
          author = line.scan(/:\s([A-Z]{2,4})_/)[0].to_s

          # Check benchmark name
          name = line.scan(/(FILECOMPARE_|BENCHMARK_)(.+)\s\.+/)[0].to_s
          name = name.gsub('FILECOMPARE_', '').gsub('BENCHMARK_', '')

          # Even test runs are benchmarks, otherwise file compares
          if (num_test_project_lines-1) % 2 == 0

            duplicate_entry = BenchmarkRun.filter(:commit_info_id => commit_info.revision, :name => name)
            if duplicate_entry.all.length > 0
              puts "Duplicate benchmark run entry"
              next
            end

            # Check benchmark time
            time = line.scan(/\s+([0-9]+\.[0-9]+)\s+sec/)[0].to_s.to_f

            benchmark_run = BenchmarkRun.create(:commit_info => commit_info,
                                                :time => time,
                                                :crashed => crashed,
                                                :name => name,
                                                :config => config,
                                                :author => Author[:short_name => author])

            #puts "Add Benchmark: #{name}, crashed #{crashed} "
          else
            # Get previous benchmark run
            bench = BenchmarkRun[:name => name, :commit_info_id => CommitInfo.order(:read_date).last.revision]
            if bench
              bench.passed = !crashed
              bench.save_changes
            else
              puts "Benchmark not found: #{line}"
            end
          end
        end
      end
    end
  end

end
