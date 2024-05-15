
if [[ -d PANTHER18.0_ascii ]];
then
    rm -rf PANTHER18.0_ascii
fi

wget ftp://ftp.pantherdb.org/panther_library/current_release/PANTHER18.0_ascii.tgz
tar -xzf PANTHER18.0_ascii.tgz
rm PANTHER18.0_ascii.tgz

