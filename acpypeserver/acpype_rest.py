import requests
from requests import get
from requests.auth import HTTPBasicAuth
import urllib.request

def AcpypePOST(user_name, password, molecule_file, charge_method,net_charge,multiplicity,atom_type):
	with requests.Session() as c:
		files = {'molecule_file': open(molecule_file,'rb')}
		url_saved = 'http://127.0.0.1:8000/submission/'
		url_login = 'http://127.0.0.1:8000/login/'
		c.get(url_login)
		csrftoken = c.cookies['csrftoken']
		login_data = dict(csrfmiddlewaretoken=csrftoken, user_name=user_name, password=password, charge_method = charge_method, net_charge = net_charge, multiplicity = multiplicity, atom_type = atom_type)
		c.post(url_login, data = dict(csrfmiddlewaretoken=csrftoken), auth=HTTPBasicAuth(user_name, password))
		r = c.post(url_saved, data= login_data, headers=dict(Referer=url_saved), files=files)
		if r.status_code == 202:
			return(print('Your submission was accepted.'))
		else:
			return(print('Your submission failed.'))

def AcpypeGETINFO(user_name, password, molecule_file):
	with requests.Session() as c:
		url_saved = 'http://127.0.0.1:8000/submission/'
		login_data = dict(user_name=user_name, password=password, molecule_file=molecule_file)
		content = c.get(url_saved, params=login_data, headers=dict(Referer=url_saved))
		data = content.text
		return data

def AcpypeGETFILES(user_name, password, func, celery_id):
	with requests.Session() as c:
		url_saved = 'http://127.0.0.1:8000/submission/'
		login_data = dict(user_name=user_name, password=password, func=func, celery_id=celery_id)
		content = c.get(url_saved, params=login_data, headers=dict(Referer=url_saved))
		data = content.text
		return data