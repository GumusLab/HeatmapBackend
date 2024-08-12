# heatviz/urls.py

from django.urls import path
from .views import process_data_view

urlpatterns = [
    path('process/', process_data_view, name='process-data'),
]
