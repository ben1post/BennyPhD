import dash_katex
import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_core_components as dcc
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import plotly
import pandas as pd

import sympy
from sympy.printing.latex import latex


def latex_out(sympyequation):
    return latex(sympyequation, order='old', long_frac_ratio = 1,mul_symbol='dot')

#'$$'+latex(sympyequation, order='old', long_frac_ratio = 1,mul_symbol='dot')+'$$'


def formula():
    """'MATHEMATIX'"""
    # MLD-dynamics
    # Material exchange between the two layers = K
    # MLD depth as function of time = M(t)
    # h+(t) is a function describing entrainment and detrainment
    # h(t) is derivtive of MLD depth
    # k is constant mixing parameter
    K, k, M, Mt, hplus, h, t = sympy.symbols('K,kappa,M,M(t),h^{+}(t),h(t),t')
    Kmix = (k+hplus)/Mt
    hplust = sympy.Max(h,0)
    ht = sympy.Derivative(Mt,t)
    K_EQ = (k+sympy.Max(ht,0))/Mt
    # For Zooplankton, K is given as K_Z = h(t) - actively maintain wihtin upper ML
    K_Z = ht/Mt

    # Nitrate
    # δD is the mineralization rate and
    # N0 is the nitrogen concentration below the upper mixed layer
    Nt,N0,deltaD_N,Dt = sympy.symbols('N(t),N_0,delta_D^N,D(t)')
    Sit,Si0 = sympy.symbols('Si(t),Si_0')
    deltaD_Si = sympy.symbols('delta_D^Si')

    NRemineralization_EQ = deltaD_N * Dt
    NMixing_EQ = K * (N0 - Nt)
    SiRemineralization_EQ = deltaD_Si * Dt
    SiMixing_EQ = K * (Si0 - Sit)

    # Phytoplankton
    # The equation describing the fitness functions of each functional type i is thus given by
    # where µP indicates the maximum growth rate and F(T) = e0.063·T is Eppley’s formulation
    # for temperature-dependent growth (Eppley, 1972).
    # The light-limiting term, H(I), represents the total light I available in
    # the upper mixed layer. According to Steele’s formulation (Steele, 1962)
    Pti, U_Ni, U_Sii, kw, PARt, OptI, one, z = \
    sympy.symbols('P_{i}(t),U^{N}_i, U^{Si}_i, k_w, PAR(t), Opt^{I}_i, 1, z')

    N_Uptake_EQ = Nt / (Nt + U_Ni)
    Si_Uptake_EQ = Sit / (Sit + U_Sii)

    LightHarvesting_EQ = (one / ((kw * Mt))) * \
        (sympy.exp(one - PARt/OptI) - (- sympy.exp(one - (PARt * sympy.exp(-kw*Mt)) / OptI)))

    I0,Iz,T = \
    sympy.symbols('I0,I_z,T')

    Iz_EQ = I0 * sympy.exp(-kw*z)
    Steele_EQ = (one/Mt * sympy.Integral(Iz/OptI * sympy.exp(1-Iz/OptI),(z,0,M)))

    TemperatureDepGrowth_EQ = sympy.exp(0.063*T)

    N_Uptake,Si_Uptake,LightHarvesting = \
    sympy.symbols('N_Uptake,Si_Uptake,LightHarvesting')

    Gains_EQ = (sympy.Min(N_Uptake_EQ, Si_Uptake_EQ)+ LightHarvesting_EQ + TemperatureDepGrowth_EQ) * Pti
    Gains = sympy.Min(N_Uptake_EQ, Si_Uptake_EQ)+ LightHarvesting + TemperatureDepGrowth_EQ

    P_Grazingi,P_Mortalityi,P_Mixingi,P_Sinkingi = \
    sympy.symbols('P^Grazing_i,P^Mortality_i,P^Mixing_i,P^Sinking_i')
    Itot,pi,R,moP,v = \
    sympy.symbols('I_tot,p_i,R,mo_P,v')

    P_Grazing = Itot * pi*Pti/R
    P_Mortality = moP
    P_Mixing = K
    P_Sinking = v / Mt

    NMixing,NRemineralization,SiMixing,SiRemineralization,Gainsi,Lossesi,i,imax = \
    sympy.symbols('N_Mixing,N_Remineralisation,Si_Mixing,Si_Remineralization,Gains_i,Losses_i,i,i_max')

    Losses_EQ = P_Grazing + P_Mortality + P_Mixing + P_Sinking
    Losses = P_Grazingi + P_Mortalityi + P_Mixingi + P_Sinkingi

    Phy = (Gainsi - Lossesi)*Pti

    Nitrate = NMixing + NRemineralization - sympy.Sum(Gainsi, (i,1,imax))
    Silicate = SiMixing + SiRemineralization - sympy.Sum(Gainsi, (i,1,imax))

    # where Is is the light level at which photosynthesis saturates and
    # I(z) is the irradiance at depth z.
    # The exponential decay of light with depth is computed according to the Beer– Lambert law
    # with a generic extinction coefficient kw

    # The current version of our model does not specify any size dependence for light absorption,
    # although we provided sug- gestions on how this could be done (Sects. 4 and 6).
    # The nutrient-limiting term U in Eq. (3) is determined by a Monod function
    # with a half-saturation constant KN, which scales allometrically
    # with phytoplankton cell size L (Litch- man et al., 2007),

    # with βU and αU, respectively, intercept and slope of the KN allometric function
    # (i.e. the power law βU ·SαU).
    # This empirical relationship is based on observations of different
    # phytoplankton groups (see Fig. 3b in Litchman et al., 2007),
    # with the regression parameters rescaled from cell volume to ESD.
    # The loss term G(Li ,Pi ) in Eq. (3) represents zooplankton grazing. As

    # The loss term V(Li ,M) in Eq. (3) represents the sinking
    # of phytoplankton as a function of size and depth of the mixed layer,

    # ZOOPLANKTON
    NMixing,deltaZ,moZ,Z_j,muZ,kz = \
    sympy.symbols('N_Mixing,delta_Z,mo_Z,Z_j,mu_Z,k_Z')

    ZooGrowth = Itot * deltaZ
    ZooMortality = moZ * Z_j ** 2
    ZooMixing = K_Z
    Zoo = ZooGrowth - ZooMortality - ZooMixing
    Itot_EQ = muZ * (R / (kz + R))
    R_EQ = sympy.Sum(pi*Pti, (i,1,imax))

    # Detritus DETRITUS DETRITUS DETRITUS DETRITUS
    # Detritus = sum(ZooMortality) + sum(UnassimilatedProduction) + sum(PhytoMortality)
    # - NRemineralization - SiRemineralization - DetritusMixing   # Detritus
    j,jmax = \
    sympy.symbols('j,j_max')

    UnassimilatedProduction = Itot * (1 - deltaZ)
    PhytoMortality = P_Mortality * Pti
    DetritusMixing_EQ = K * Dt

    Detritus = sympy.Sum(ZooMortality+UnassimilatedProduction,(j,1,jmax)) + \
    sympy.Sum(PhytoMortality,(i,1,imax)) - NRemineralization_EQ - SiRemineralization_EQ - DetritusMixing_EQ

    # initialize List of List of all relevant formula and parameters
    formulaframe = [['Nitrate', latex_out(Nt), latex_out(Nitrate)],
                    ['Nitrate', latex_out(Nt), latex_out(Nitrate)],
                    ['Nitrate', latex_out(Nt), latex_out(Nitrate)]]

    # Create the pandas DataFrame
    df = pd.DataFrame(formulaframe, columns=['State Variable', 'Symbol','Equation'])

    return df

#' TODO:
#' so #1 I need all the equations, as simply the process names, and additionally the actual mathmatical formula!
#' this should be stored inside a pandas dataframe, so that I can be dynamically calles by the dash app
#' now I have to start implementing that!
#'

dfmath = formula()

s = dfmath[dfmath['State Variable'] == 'Nitrate']#hoverData['points'][0]['customdata']]

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
    dcc.Graph(
        id='ODE',
        figure={
            'data': [{
                    'x': [1, 2, 3],
                    'y': [1, 2, 3],
                    'Name': dfmath['State Variable'],
                    'text': dash_katex.DashKatex(
                            id = 'eqtest',
                            expression='The {}, {} {} '.format(
                                s.iloc[0]['State Variable'],
                                s.iloc[0]['Symbol'],
                                s.iloc[0]['Equation'])),
                    'mode': 'text+EQ',
                    'hovertext':dfmath['State Variable'],
            'marker': {
                'color': 'grey',
                'size': 8,
                'opacity': 0.6
            },
            'customdata': ['Equation'],
            #'type': 'scattermapbox',
            }],
            'layout': {
                'title': 'Dash Data Visualization',
                'hovermode': 'closest',
                'margin': {'l': 0, 'r': 0, 'b': 0, 't': 0},  }
            }),

    html.Div('Here is some latex:'),
    dash_katex.DashKatex(
        id='katex',
        expression=""
        )
    ])



@app.callback(Output('katex', 'expression'),
              [Input('ODE', 'hoverData')])
def updateExpression(hoverData):
    return hoverData


@app.callback(Output('katex', 'displayMode'),
              [Input('radio-item-display-mode', 'value')])
def updateDisplayMode(value):
    return value


@app.callback(Output('katex', 'throwOnError'),
              [Input('radio-item-throw-on-error', 'value')])
def updatethrowOnError(value):
    return value


if __name__ == '__main__':
    app.run_server(debug=True)
