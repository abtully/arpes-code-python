from dash_dashboard.app import get_app
from ali_dash import leed_analysis_page

app = get_app([leed_analysis_page])

if __name__ == '__main__':
    app.run_server(port=8060, debug=True)