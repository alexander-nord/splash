
if [[ -d PANTHER18.0 ]];
then
    rm -rf PANTHER18.0
fi

wget ftp://ftp.pantherdb.org/panther_library/current_release/PANTHER18.0_ascii.tgz
tar -xzf PANTHER18.0_ascii.tgz

mv target/famlib/rel/PANTHER18.0_altVersion/ascii/PANTHER18.0 .

rm PANTHER18.0_ascii.tgz

