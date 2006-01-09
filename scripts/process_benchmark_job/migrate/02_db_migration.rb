# To run that execute 'sequel -m . sqlite:///home/lars/db/ogsbench_test.db'

Sequel.migration do
  up do
    add_column :commit_infos, :is_svn_commit, Fixnum
    add_column :commit_infos, :read_date, Time
    self[:commit_infos].update(:is_svn_commit=>'1')
    self[:commit_infos].update(:read_date=>Time.now)
  end
  down do
    drop_column :commit_infos, :is_svn_commit
    drop_column :commit_infos, :read_date
  end
end