from django.urls import path
from .views import process_data_view, cleanup_session,command_execution,correlation_network,refresh_heatmap,estimate_correlation_job,recommend_correlation_parameters,load_example_data,enrich_analysis_view,get_3d_coords_view

urlpatterns = [
    path('process/', process_data_view, name='process-data'),
    path('cleanup/', cleanup_session, name='cleanup-session'),  
    path('command/', command_execution, name='command-execution'),
    path('network-correlation/', correlation_network, name='correlation-network'),
    path('api/estimate-correlation-job/', estimate_correlation_job, name='estimate_correlation_job'),
    path('api/recommend-correlation-parameters/', recommend_correlation_parameters, name='recommend_correlation_parameters'),
    path('refresh-heatmap/', refresh_heatmap, name='refresh-heatmap'),
    path('examples/load/', load_example_data, name='load-example-data'), 
    path('enrich/', enrich_analysis_view, name='enrich_analysis'),

     # NEW: URL for fetching 3D coordinates
    path('get_3d_coords/<str:session_id>/', get_3d_coords_view, name='get_3d_coords'),
]
