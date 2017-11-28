import django_tables2 as tables
from .models import Submition

class SimpleTable(tables.Table):
    class Meta:
        model = Submition