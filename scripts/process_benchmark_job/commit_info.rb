require 'rubygems'
require 'time'
require 'sequel'
#require 'database.rb'
require './ogs_author_mapping.rb'

Sequel::Model.unrestrict_primary_key

# Create a table if it not exists
$DB.create_table? :commit_infos do
  primary_key     :revision
  String          :git_revision
  Fixnum          :is_svn_commit
  Time            :date
  Time            :read_date
  String          :branch
  foreign_key     :author_id, :table => :authors
end

class CommitInfo < Sequel::Model(:commit_infos)
  set_primary_key :revision
  set_dataset dataset.order(:revision)
  many_to_one :author
  many_to_one :benchmark_run
end

class CommitInfoLoader

  def new?
    @new
  end
  
  attr_reader :commit_info

  def initialize(filename)

    @new = true

    puts filename
    File.open(filename, 'r') do |file|
      
      svn = false
      git = false
      first_line = file.gets
      if first_line =~ /Path:/
        svn = true
      elsif first_line =~ /commit/
        git = true
      end
      
      revision = 0
      git_revision = ""
      author = nil
      date = nil
      read_date = Time.now
      branch = nil
      is_svn_commit = 0

      if svn
        is_svn_commit = 1
        while line = file.gets
          line.scan(/svn\/ogs\/([\S]+)\//) do |match|
            branch = match[0]
          end
          line.scan(/Revision:\s([0-9]+)/) do |match|
            revision = match[0].to_i
          end
          line.scan(/Last Changed Author:\s([\w]+)/) do |match|
            author_name = match[0]
            author = Author[:svn_user => author_name]
          end
          line.scan(/Last Changed Date:\s([0-9]{4}-[0-9]{2}-[0-9]{2}\s[0-9]{2}:[0-9]{2}:[0-9]{2})/) do |match|
            date = Time.parse(match[0])
          end
        end
      elsif git
        first_line.scan(/commit ([\S]{40,40})/) do |match|
          git_revision = match[0]
          # Get only first 12 chars and convert to integer
          revision = git_revision[0,12].to_i(16)
        end
        while line = file.gets
          line.scan(/Author:\s([\S\s]+)\s</) do |match|
            author = Author[:name => match[0]]
          end
          line.scan(/Date:[\s]+([\s\S]+)/) do |match|
            date = Time.parse(match[0])
          end
        end
      end

      if CommitInfo[:revision => revision]
        @new = false
	      if svn
          puts "Commit info of revision #{revision} already read."  
        else
	        puts "Commit info of git revision #{git_revision} and revision #{revision} already read."
	      end
      else
        puts revision
        puts git_revision
        puts read_date
        @commit_info = CommitInfo.create(:revision => revision,
                                        :git_revision => git_revision,
                                        :date => date,
                                        :read_date => read_date,
                                        :branch => branch,
                                        :is_svn_commit => is_svn_commit)
        @commit_info.author = author
        @commit_info.save
      end

    end
  end

end

#CommitInfoLoader.new('tests/svnInfoOld.txt')
#$DB[:commit_infos].each {|row| p row}
