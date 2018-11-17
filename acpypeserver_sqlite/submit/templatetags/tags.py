from django import template
from django.conf import settings

register = template.Library()


@register.filter
def statusColor(value):
    color = {'Deleted':'orange', 'Queued':'black','Running':'blue', 'Failed':'red',
              'Cancelled':'magenta', 'Finished':'green'}
    return color.get(value, 'black')


@register.simple_tag
def settings_value(name):
    return getattr(settings, name, "")