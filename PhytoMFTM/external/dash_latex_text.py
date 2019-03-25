import Tests.
import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_core_components as dcc

app = dash.Dash(__name__)

app.scripts.config.serve_locally = True
app.css.config.serve_locally = True

app.layout = html.Div([
    html.Div("Expression Input:"),
    dcc.Input(id='input-latex-expression'),

    html.Div("Display Mode:"),
    dcc.RadioItems(
        id='radio-item-display-mode',
        options=[
            {'label': 'True (Display/Equation Mode)', 'value': True},
            {'label': "False (Inline Mode)", 'value': False}
        ]
    ),

    html.Div("Throw on error:"),
    dcc.RadioItems(
        id='radio-item-throw-on-error',
        options=[
            {'label': 'True', 'value': True},
            {'label': "False", 'value': False}
        ]
    ),

    html.Div('Here is some latex:'),
    dash_katex.DashKatex(
        id='katex',
        expression=""
    ),
    html.Div(id='main',children=[
    html.Div(id='static',children='$$ x=1 $$'),
    html.Div(id='dynamic',children=''),
    html.Button('Add Math',id='button'),
])
])


@app.callback(Output('katex', 'expression'),
              [Input('input-latex-expression', 'value')])
def updateExpression(value):
    return value


@app.callback(Output('katex', 'displayMode'),
              [Input('radio-item-display-mode', 'value')])
def updateDisplayMode(value):
    return value


@app.callback(Output('katex', 'throwOnError'),
              [Input('radio-item-throw-on-error', 'value')])
def updatethrowOnError(value):
    return value



app = dash.Dash(__name__)
#mathjax = '/Users/bp/Desktop/ToyModels/BennyPhD/PhytoMFTM/Tests/MathJax.js'#
#app.scripts.append_script({ mathjax })


@app.callback(
  dash.dependencies.Output('dynamic','children'),
 [dash.dependencies.Input('button','n_clicks')]
)
def addmath(n_clicks):
  if n_clicks:
    return '$$ x=1 $$'
  else:
    return None

if __name__ == '__main__':
  app.run_server(debug=True)