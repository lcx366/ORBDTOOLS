from astropy.time import Time
from odtools.transform.frame_trans import gcrf_teme_mat

stations_url = 'TLE_p0/tle_20220524.txt'
satellites = load_tle_file(stations_url)
times = Time(['2023-07-01T09:07:38','2023-06-27T16:35:22', '2023-06-29T04:35:22','2023-06-26T13:35:22','2023-06-28T04:35:22','2023-07-01T09:07:38','2023-06-27T16:35:22', '2023-06-29T04:35:22','2023-06-26T13:35:22','2023-06-28T04:35:22'])

import time

t_start = time.time()

e, r_teme, v_teme = satellites.sgp4(times.jd1, times.jd2)
gcrf2teme_mat,teme2gcrf_mat = gcrf_teme_mat(times)
r_gcrf = (teme2gcrf_mat @ r_teme.transpose(1,2,0)).transpose(2,0,1)

t_end = time.time()
print('time used:',t_end-t_start)

from odtools.targetmatch.target_matching import targetid_find

import time

t_start = time.time()
status = targetid_find(out_filename,tle_file,mode,logger=satmatch_logger)
print(status)
print('time used:',t_end-t_start)
