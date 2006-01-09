require 'rubygems'
require 'sequel'
#require 'database.rb'


# Create a table if it not exists
$DB.create_table? :authors do
  primary_key :id
  String :name
  String :svn_user
  String :email
  String :short_name
  unique :svn_user
end


class Author < Sequel::Model(:authors)
  one_to_one :commit_info
  one_to_one :benchmark_run
end

class AuthorLoader
  
  def initialize(file_name)
    File.open(file_name, 'r') do |file|

      while line = file.gets
        # Ignore commented lines
        if line =~ /^#/
          next
        end

        svn_user = line.scan(/^[\w]+/).to_s
        name = line.scan(/=.+</).to_s
        name = name.gsub(/^=\s/, '').gsub(/\s</, '')
        email = line.scan(/[A-Za-z0-9+_.-]+@[A-Za-z0-9.-]+/).to_s
        short_name_match = name.gsub(/Prof.\s|Dr.\s/, '').scan(/[A-Z]|\s[a-z]/)
        short_name = ''
        short_name_match.each { |name_char|
          short_name = short_name + name_char.gsub(/\s/, '').upcase
        }

        if Author[:svn_user => svn_user]
          puts "Author #{Author[:svn_user => svn_user].name} already registered."
        else
          # Insert into database
          Author.create(:name => name,   :svn_user => svn_user,
                        :email => email, :short_name => short_name)
        end
      end
    end
  end
end


#AuthorLoader.new('authors.txt')
#$DB[:authors].each {|row| p row}
#print $DB[:authors].length

