from django import template

register = template.Library()

@register.filter
def splitx(value, arg=" "):
    lista = value.split(arg)
    lista = map(str,lista)
    return lista

@register.filter
def pos(value, arg = 0):
    return value[int(arg)]

@register.filter
def bname(value):
    import os
    return os.path.basename(value)

@register.filter
def statusColor(value):
    color = {'Submitted':'orange', 'Running':'blue', 'Failed':'red',
              'Cancelled':'magenta', 'Finished':'green'}
    return color.get(value, 'black')

@register.filter
def replace(value, arg):
    arg1, arg2 = arg.split()
    return value.replace(arg1,arg2)

@register.filter
def objdate(value):
    import datetime
    return datetime.datetime.strptime(str(value), '%a %b %d %H:%M:%S %Y')
