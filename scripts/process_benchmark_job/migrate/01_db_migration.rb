# To run that execute 'sequel -m . sqlite:///home/lars/db/ogsbench_test.db'

Sequel.migration do
  up do
    add_column :commit_infos, :git_revision, String
    self[:commit_infos].update(:git_revision=>'')
  end
  down do
    drop_column :commit_infos, :git_revision
  end
end