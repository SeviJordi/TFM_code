# Jordi Sevilla Fortuny
# Script to complete the install non conda dependecies

cd $HOME

sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

echo "export PATH=\$HOME/edirect:\$PATH" >> $HOME/.bash_profile

exit 0