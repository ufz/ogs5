#!/usr/bin/env ruby

require 'rubygems'
require 'sequel'

$debug = true

if ARGV.size > 4
  $DB = Sequel.connect('sqlite://' + ARGV[4])
  puts "Connecting to #{ARGV[4]}"
else
  $DB = Sequel.connect('sqlite://'+ ENV['HOME'] + '/db/ogsbench.db')
end

require './bench_info.rb'
require './commit_info.rb'
require './ogs_author_mapping.rb'
require 'net/smtp'
require 'csv'

class BenchmarkInfoProcessor

  # Use BenchmarkInfoProcessor(new_revision) or BenchmarkInfoProcessor(new_revision, old_revision)
  def initialize(commit_info)

    @actual_benchmark_runs = []
    @last_benchmark_runs = []
    @failed_benchmarks = []
    @crashed_benchmarks = []
    @new_failed_benchmarks = []
    @fixed_benchmarks = []

    @new_commit_info = commit_info
    if commit_info.is_svn_commit == 1
      @old_commit_info = CommitInfo.filter('revision < ?', @new_commit_info.revision).order(:revision).last
    else
      @old_commit_info = CommitInfo.filter('read_date < ?', @new_commit_info.read_date).order(:read_date).last
    end

    puts "Comparing benchmarks of revision #{@old_commit_info.revision} and revision #{@new_commit_info.revision}:" if $debug
    process
  end


  def process

    ## Get failed benchmarks
    @actual_benchmark_runs = BenchmarkRun.filter('commit_info_id => ?', @new_commit_info.revision)
    @failed_benchmarks = @actual_benchmark_runs.filter(:passed => false)
    @crashed_benchmarks = @failed_benchmarks.filter(:crashed => true)

    ## Get newly failed benchmarks
    # Get last benchmark run
    @last_benchmark_runs = BenchmarkRun.filter(:commit_info_id => @old_commit_info.revision)


    @failed_benchmarks.each do |actual_failed_benchmark|
      # If it is completely new then add
      if @last_benchmark_runs.filter(:name => actual_failed_benchmark.name).all.length == 0
        @new_failed_benchmarks.push actual_failed_benchmark
      # Otherwise check if it ran successfully
      else
        new_failed_benchmark = @last_benchmark_runs.filter(:passed => true,
                                                          :name => actual_failed_benchmark.name).all
        if new_failed_benchmark.length > 0
          @new_failed_benchmarks.push new_failed_benchmark[0]
      end
      end
    end

    ## Get fixed benchmarks

    @last_benchmark_runs.filter(:passed => false).each do |last_benchmark|
      #puts "last failed: #{last_benchmark.name}"
      fixed_benchmark = @actual_benchmark_runs.filter(:passed => true,
                                                      :name => last_benchmark.name).all
      if fixed_benchmark.length > 0
        @fixed_benchmarks.push fixed_benchmark[0]
      end
    end

  end

  def print_summary()
    ## Generate statistics
    num_tested = @actual_benchmark_runs.count
    num_failed = @failed_benchmarks.count

    puts ""
    puts "########## BENCHMARK RUN SUMMARY ##########"
    if @fixed_benchmarks.length == 0
      puts "No fixed benchmarks."
    else
      puts "Fixed benchmarks:"
      @fixed_benchmarks.each {|row| p row.name}
    end

    if @new_failed_benchmarks.length == 0
      puts "No new failed benchmarks."
    else
      puts "New failed benchmarks:"
      @new_failed_benchmarks.each {|row| p row.name}
    end

    puts "Summary of all run benchmarks: #{((num_tested.to_f - num_failed) / num_tested * 100).to_i} % passed"
    puts "###########################################"
    puts ""
  end

  def write_statistics_to_csv(filename)
    ## Generate statistics
    num_tested = @actual_benchmark_runs.count
    num_failed = @failed_benchmarks.count
    num_new_failed = @new_failed_benchmarks.length
    num_crashed = @crashed_benchmarks.count
    num_fixed = @fixed_benchmarks.length

    # Write to csv file
    #CSV.open("#{File.dirname(ARGV[0])}/benchSummary.csv", "w") do |csv|
    CSV.open(filename, "w") do |csv|
      csv << ['tested', 'failed', 'crashed', 'new_failed', 'fixed']
      csv << [num_tested.to_s, num_failed.to_s, num_crashed.to_s,
              num_new_failed.to_s, num_fixed.to_s]
    end
  end

  def send_email()
    BenchmarkEmail.new(@new_failed_benchmarks, @fixed_benchmarks, @failed_benchmarks)
  end

end

class BenchmarkEmail

  def initialize(new_failed_benchmarks, fixed_benchmarks, failed_benchmarks)

    # Write email to commiter if something has changed
    if new_failed_benchmarks.length > 0 or fixed_benchmarks.length > 0
      author = Author[:id => CommitInfo.order(:read_date).last.author.id]
      msg = "Hello #{author.name},\n\n"

      nice_verbs = ['awesome', 'brilliant', 'great', 'gorgeous']

      if fixed_benchmarks.length > 0
        msg << "you are absolutely #{nice_verbs[rand(nice_verbs.length)]}! You have fixed the following benchmarks:\n"
        fixed_benchmarks.each do |benchmark|
          msg << benchmark.benchmark_info << "\n"
        end
        msg << "\n"
      end

      bla = ['Ah', 'Eh', 'Ouch']

      if new_failed_benchmarks.length > 0
        msg << "#{bla[rand(bla.length)]}, something went wrong because there are new benchmarks failing:\n"
        new_failed_benchmarks.each do |benchmark|
          msg << benchmark.benchmark_info << "\n"
          msg << "Have a look at #{$job_url.to_s + benchmark.author.short_name + "_" + benchmark.name.gsub(/\//, "_") + ".html"} for more infos\n\n"
        end
        msg << "\nPlease fix them again as soon as possible.\n"

      end

      if failed_benchmarks.count > 0
        msg << "Unfortunately the following benchmarks failed as before:\n"
        failed_benchmarks.each do |benchmark|
          if not new_failed_benchmarks.find{ |failed_benchmark| failed_benchmark.name == benchmark.name}
            msg << benchmark.benchmark_info << "\n"
            msg << "Have a look at #{$job_url.to_s + "artifact/benchmarks/results/" + benchmark.author.short_name + "_" + benchmark.name.gsub(/\//, "_") + ".html"} for more infos\n\n"
          end
        end
      end

      puts "Sending email to #{author.name} <#{author.email}>"
      #send_email('Hudson Build Server', 'lars.bilke@ufz.de',
      #           author.email, author.name,
      #           'Benchmark report', msg) if $password
      send_email('Hudson Build Server', 'lars.bilke@ufz.de',
                 'lars.bilke@ufz.de', author.name,
                 'Benchmark report', msg) if $password
    end

  end

  # if using text/html as Content-type then newlines must be encoded as <br> ...
  def send_email(from, from_alias, to, to_alias, subject, message)
  	msg = <<END_OF_MESSAGE
From: #{from_alias} <#{from}>
To: #{to_alias} <#{to}>
MIME-Version: 1.0
Content-type: text/plain
Subject: #{subject}

#{message}
END_OF_MESSAGE

  	Net::SMTP.start('imap.leipzig.ufz.de',
                    25,
                    'localhost.localdomain',
                    'bilke', $password, :plain) do |smtp|
  		puts "Sending email to #{to_alias} <#{to}>..."
      smtp.send_message(msg, from, to)
  	end
  end

end

### Main ###
if false
  bi = BenchmarkInfoProcessor.new
  #bi.write_statistics_to_csv('test.csv')
  #bi.send_email
  bi.print_summary
else

  if ARGV.size < 2
    puts 'Usage: process_benchmark_job.rb commit_info_file benchmark_job_output'
    Process.exit 1
  end

  $password = nil
  $job_url = nil

  if ARGV.size > 2
    $password = ARGV[2]
  end

  if ARGV.size > 3
    $job_url = ARGV[3]
  end

  ### Read files ###
  ci = nil
  if File.exists?(ARGV[0])
    # read commit info
    ci = CommitInfoLoader.new(ARGV[0])

    if File.exists?(ARGV[1])
      # read benchmark job output
      if  ci.commit_info
        BenchmarkRunsLoader.new(ARGV[1], ci.commit_info)
      else
        puts "Aborting: Commit already processed."
        Process.exit 1
      end
    else
      puts "File #{ARGV[1]} does not exist!"
      Process.exit 1
    end
  else
    puts "File #{ARGV[0]} does not exist!"
    Process.exit 1
  end

  ### Process info ###
  bi = BenchmarkInfoProcessor.new(ci.commit_info)
  bi.write_statistics_to_csv("#{File.dirname(ARGV[0])}/benchSummary.csv")
  bi.send_email if ci.new?
  bi.print_summary

  ## Generate plots ##
  if false
    bench_dirs = [ 'dir1/', 'dir2/' ]
    benchmark_dir = './../../benchmarks/'
    reference_dir = './../../benchmarks_ref/'
    pdflatex_dir = '/usr/bin'
    tioga_dir = '/usr/bin'

    # Cleanup from previous runs
    Dir.glob("./*.pdf") do |file|
      File.delete(file)
    end

    bi.failed_benchmarks.each do |actual_failed_benchmark|
      bench_dir = File.dirname(actual_failed_benchmark.name) + "/"
      puts "Searching for tec files in " + benchmark_dir + bench_dir
      next if not File.directory?(benchmark_dir + bench_dir)
      # Search for .tec files
      Dir.glob(benchmark_dir + bench_dir + '*.tec') do |new_tecfile|
        #puts new_tecfile
        # Check for reference tecfile
        if File.exists?(new_tecfile.gsub(benchmark_dir, reference_dir))
          ref_tecfile = new_tecfile.gsub(benchmark_dir, reference_dir)

          # Create input file for tioga
          tioga_input_file = File.read('tioga_main_template.rb')
          tioga_input_file.gsub!('$ref_tecfile$', ref_tecfile)
          tioga_input_file.gsub!('$new_tecfile$', new_tecfile)

          output_filename = new_tecfile.gsub(benchmark_dir, "").gsub(/\//, "_")
          tioga_input_filename = "#{output_filename}.rb"
          File.open(tioga_input_filename, "w") do |file|
            file.puts tioga_input_file
          end

          # Run tioga
          %x[export PATH=#{pdflatex_dir}:$PATH && #{tioga_dir}/tioga #{tioga_input_filename} -p]

          # Delete tioga input file
          #File.delete(tioga_input_filename)

          portfolio_filename = output_filename + "_portfolio.pdf"
          if File.exists?(portfolio_filename)
            # Rename plot files
            File.rename(portfolio_filename, portfolio_filename.gsub("_portfolio", ""))

            # Output
            puts "SUCCESS: Plot written to #{portfolio_filename.gsub("_portfolio", "")}"
          else
            puts "ERROR: Something went wrong when running tioga! No plot file was written."
          end

        else
          # If it not exists goto next tecfile
          puts "ERROR: Reference file #{new_tecfile.gsub(benchmark_dir, reference_dir)} missing!"
          next
        end


      end
    end

    # Delete intermediate plot files ending in plot.pdf
    Dir.glob("./*plot.pdf") do |file|
      File.delete(file)
    end

  end # Generate plots

end # if true

$DB.disconnect
