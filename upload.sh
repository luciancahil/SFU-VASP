

username=royhe62
main_remote=/home/royhe62/projects/def-samiras-ab/royhe62/SFU-VASP

scp -r processing/*.traj ${username}@cedar.alliancecan.ca:${main_remote}/processing

scp -r POSCARS/*.csv ${username}@cedar.alliancecan.ca:${main_remote}/POSCARS
