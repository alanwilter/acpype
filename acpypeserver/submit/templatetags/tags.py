from django import template

register = template.Library()

@register.filter
def statusColor(value):
    color = {'Deleted':'orange', 'Running':'blue', 'Failed':'red',
              'Cancelled':'magenta', 'Finished':'green'}
    return color.get(value, 'black')
