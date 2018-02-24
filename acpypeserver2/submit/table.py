import django_tables2 as tables
from .models import Submission

class SimpleTable(tables.Table):
    class Meta:
        model = Submission