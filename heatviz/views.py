# heatviz/views.py

from django.shortcuts import render
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status
from rest_framework.parsers import MultiPartParser
from rest_framework.decorators import parser_classes
from heatviz.clustergrammer2.make_clustergrammer import make_cluster
import pandas as pd

@api_view(['POST'])
@parser_classes([MultiPartParser])
def process_data_view(request):

    file = request.FILES.get('data')  # Get the file from the request
    if not file:
        return Response({"error": "No file provided"}, status=status.HTTP_400_BAD_REQUEST)
    try:
        data = pd.read_csv(file,index_col=0)  # Read the file using pandas
        response_data = make_cluster(data)
        return Response(response_data, status=status.HTTP_200_OK)
    except Exception as e:
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
