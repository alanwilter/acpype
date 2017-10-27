from submit.models import Submition

def run():
	all_products=Submition.objects.all()
	print all_products