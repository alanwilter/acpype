from submit.models import Submission
from rest_framework import serializers


class SubmissionSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Submission
        fields = ('molecule_file', 'charge_method', 'net_charge', 'multiplicity', 'atom_type', 'juser', 'jstatus', 'jzipped', 'jlog', 'date', 'usr_folder')