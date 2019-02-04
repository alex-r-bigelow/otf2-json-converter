phylanx-data-wrangling
======================
Some wrangling scripts for OTF2 + tree data from running phylanx. Ideally, [Origraph](https://origraph.github.io) should be able to handle this, but:
1. There are non-graph wrangling operations involved, and
2. The OTF2 traces are way too big to fit in memory (and need special parsing)

Setup
=====
```bash
git clone https://github.com/alex-r-bigelow/phylanx-data-wrangling.git
cd phylanx-data-wrangling
sudo apt install python3-venv
python3 -m venv env
source env/bin/activate
pip3 install -r requirements.txt
```
(TODO: generate requirements.txt)

Example
=======
```bash
./convertToJson.py --input test_data/stdout.txt --events >results.json
```