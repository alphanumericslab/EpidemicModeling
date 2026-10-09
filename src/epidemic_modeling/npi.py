"""Matched MATLAB/Python training, forecasting, and NPI control.

This explicit, deterministic API replaces the original plotting-heavy
XPRIZE orchestration. The SI-alpha and optimal-control equations are retained.
Author: Reza Sameni, Emory University. doi:10.1109/JSTSP.2021.3129118.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .data import prepare_cases, read_oxford_data
from .kalman import si_alpha_model_ekf, si_alpha_model_ekf_opt_controlled
from .models import si_alpha_controlled

NPI_COLUMNS = ('C1_School closing','C2_Workplace closing','C3_Cancel public events',
    'C4_Restrictions on gatherings','C5_Close public transport','C6_Stay at home requirements',
    'C7_Restrictions on internal movement','C8_International travel controls',
    'H1_Public information campaigns','H2_Testing policy','H3_Contact tracing','H6_Facial Coverings')
NPI_MAXES = np.array([3,3,2,4,2,3,2,4,2,3,2,4],float)


def npi_cost(new_cases, inputs, weights):
    """Return mean daily cases and mean elementwise weighted NPI intensity.

    inputs is interventions-by-time; weights is a scalar, intervention
    vector, or matching matrix. The denominator for NPI cost includes both
    interventions and days, exactly as in the original MATLAB function.
    """
    u = np.asarray(inputs,float); w = np.asarray(weights,float)
    if w.ndim == 1: w = w[:,None]
    return float(np.mean(new_cases)),float(np.mean(w*u))


def nonnegative_least_squares(x,y,max_iter=10000,tolerance=1e-10):
    """Solve nonnegative least squares by matched cyclic coordinate descent.

    Both languages use the same zero initialization, coordinate order,
    relative stopping criterion, and iteration cap. This avoids differences
    between MATLAB's lsqnonneg and SciPy solvers in the updated pipeline.
    """
    x = np.asarray(x,float); y = np.asarray(y,float).reshape(-1)
    if x.ndim != 2 or x.shape[0] != y.size or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("finite x/y with matching rows required")
    coef = np.zeros(x.shape[1]); residual = y.copy(); norm = np.sum(x*x,axis=0)
    for _ in range(max_iter):
        old = coef.copy()
        for j in range(len(coef)):
            if norm[j] == 0: continue
            delta = max(0,coef[j]+x[:,j]@residual/norm[j])-coef[j]
            coef[j] += delta; residual -= delta*x[:,j]
        if np.max(np.abs(coef-old),initial=0) <= tolerance*(1+np.max(np.abs(coef),initial=0)): break
    return coef


def default_si_params(u_max, *, epsilon=.5, weights=None):
    """Create reproducible SI-alpha parameters used by both implementations.

    beta=-log(.01)/21, gamma=1/7, dt=1 day; other choices are documented
    educational assumptions from the historical implementation.
    """
    maximum = np.asarray(u_max,float).reshape(-1)
    return dict(dt=1.,beta=-np.log(.01)/21,gamma=1/7,a=np.zeros(len(maximum)),b=0.,
        u_min=np.zeros(len(maximum)),u_max=maximum,alpha_min=0.,alpha_max=5.,
        epsilon=epsilon,w=np.ones(len(maximum)) if weights is None else np.asarray(weights,float),
        sigma=10000.,obs_type='NEWCASES',s_min=0.,i_min=0.)


def fit_npi_model(cumulative, inputs, population, u_max, *, regression_start=0):
    """Fit a two-pass EKF/NNLS SI-alpha model to one geographic time series.

    cumulative is length N; inputs is (interventions,N); population is a
    positive count. The first pass estimates alpha with zero NPI response;
    NNLS regresses alpha on u_max-u, and the second pass refines alpha with
    that response. A second NNLS uses the refined series (fixing the legacy
    first-series reuse). Returns a serializable dict with parameters, final
    state/covariance, coefficients and training length. See docs/migration.md.
    """
    cumulative = np.asarray(cumulative,float).reshape(-1); u = np.asarray(inputs,float)
    maximum = np.asarray(u_max,float).reshape(-1)
    if len(cumulative)<14 or population<=0 or u.shape!=(len(maximum),len(cumulative)):
        raise ValueError("need >=14 days, population>0, and aligned interventions")
    if not 0<=regression_start<len(cumulative): raise ValueError("invalid regression_start")
    if not np.all(np.isfinite(u)) or np.any(u<0) or np.any(u>maximum[:,None]):
        raise ValueError("interventions must be finite and within bounds")
    daily,smoothed = prepare_cases(cumulative)
    positives = smoothed[smoothed>0][:7]; initial = max(10,float(positives.mean()) if positives.size else 10)
    if initial>=population: raise ValueError("population must exceed initial case count")
    fraction = initial/population; params = default_si_params(maximum)
    state = np.array([1-fraction,fraction,params['beta']+np.log(2.5)])
    covariance = 100*np.diag([fraction,fraction,.01])**2
    q = np.diag([10*fraction,30*fraction,.01])**2
    variance = max(float(np.var((daily-smoothed)/population,ddof=1)),1e-14)
    common = [state,covariance,np.full(3,np.nan),np.full((3,3),np.nan),np.zeros(3),0.,q,variance,1.,1.,21,1]
    first = si_alpha_model_ekf(np.zeros_like(u),smoothed/population,params,*common)
    design = (maximum[:,None]-u).T
    coef = nonnegative_least_squares(design[regression_start:],first.s_smooth[2,regression_start:])
    params['a'] = coef
    second = si_alpha_model_ekf(u,smoothed/population,params,*common)
    refined = nonnegative_least_squares(design[regression_start:],second.s_smooth[2,regression_start:])
    params['a'] = refined
    return dict(population=float(population),params={k: v.tolist() if isinstance(v,np.ndarray) else v for k,v in params.items()},
        coefficients=coef.tolist(),refined_coefficients=refined.tolist(),state=second.s_plus[:,-1].tolist(),
        covariance=second.p_plus[:,:,-1].tolist(),training_days=len(cumulative))


def forecast_npi(model, inputs):
    """Forecast daily cases and post-update states under a supplied NPI schedule.

    Returns (cases, states), where states is (3,N). Cases use the pre-update
    incidence N_population*s*i*alpha, avoiding a one-day shift between
    forecasted cases and the transition producing those cases.
    """
    u = np.asarray(inputs,float); p = model['params']; initial = np.asarray(model['state'],float)
    if u.ndim == 1: u = u.reshape(1,-1)
    if u.shape[0] != len(p['a']) or not np.all(np.isfinite(u)) or np.any(u<np.asarray(p['u_min'])[:,None]) or np.any(u>np.asarray(p['u_max'])[:,None]):
        raise ValueError("forecast controls must be finite and within model bounds")
    states = np.array(si_alpha_controlled(u,*initial,p['u_max'],p['alpha_min'],p['alpha_max'],
        p['gamma'],p['a'],p['b'],p['beta'],0,0,0,u.shape[1],p['dt'],noise=np.zeros((3,u.shape[1]))))
    before = np.column_stack((initial,states[:,:-1]))
    return model['population']*np.prod(before,axis=0),states


def optimal_npi(model, horizon, weights, epsilon=.5, iterations=8):
    """Solve a bounded finite-horizon SI-alpha control problem by EKF/EKS shooting.

    Six states include three costates. Repeated filtering/smoothing imposes
    zero terminal costates; NaN controls select lower/upper bounds using
    the paper's switching function. Returns (controls, cases, states).
    This is a finite-iteration solver, not a claim of globally optimal convergence.
    """
    if not 0<=epsilon<=1 or horizon<1 or int(horizon)!=horizon or iterations<1:
        raise ValueError("epsilon in [0,1], integer horizon>=1, iterations>=1 required")
    p = dict(model['params']); p.update(epsilon=epsilon,w=np.asarray(weights,float),sigma=10000.)
    for key in ('a','u_min','u_max'): p[key] = np.asarray(p[key],float)
    if p['w'].shape != p['a'].shape or np.any(p['w']<0): raise ValueError("nonnegative weights must match interventions")
    initial = np.r_[model['state'],np.zeros(3)]; cov = np.zeros((6,6))
    cov[:3,:3] = model['covariance']; cov[3:,3:] = np.eye(3)
    terminal = np.r_[np.full(3,np.nan),np.zeros(3)]
    terminal_cov = np.full((6,6),np.nan); terminal_cov[3:,3:] = 0
    u = np.full((len(p['a']),horizon),np.nan)
    for _ in range(iterations):
        result = si_alpha_model_ekf_opt_controlled(u,np.full(horizon,np.nan),p,initial,cov,terminal,terminal_cov,
            np.zeros(6),0,np.diag([1e-10,1e-10,1e-6,1e-4,1e-4,1e-4]),1e-6,1,1,21,1)
        initial[3:] = result.s_smooth[3:,0]
    controls = result.u_opt_smooth.copy(); controls[:,-1] = result.u_opt[:,-1]
    cases,states = forecast_npi(model,controls)
    return controls,cases,states


def train_npi_prescriptor(start_date_str,end_date_str,data_file,geo_file,populations_file,included_ip,npi_maxes,trained_model_params_file, *, regression_start_date=None):
    """Train selected geographies from local Oxford data and save a JSON model bundle.

    CountryName/RegionName identify rows in all three input CSVs.
    Population2020 is required in the population table. Dates are inclusive
    ISO strings. Missing interventions forward-fill then become zero.
    Returns the same bundle written to trained_model_params_file.
    """
    data = read_oxford_data(data_file); geos = pd.read_csv(geo_file).fillna({'RegionName':''})
    populations = pd.read_csv(populations_file).fillna({'RegionName':''}); models = []
    for geo in geos.itertuples(index=False):
        country,region = geo.CountryName,geo.RegionName
        segment = data[(data.CountryName==country)&(data.RegionName==region)&
            (data.Date>=pd.Timestamp(start_date_str))&(data.Date<=pd.Timestamp(end_date_str))]
        if len(segment)<14: continue
        if not (segment.Date.diff().dropna()==pd.Timedelta(days=1)).all():
            raise ValueError(f"daily observations must be consecutive: {country}/{region}")
        pop = populations[(populations.CountryName==country)&(populations.RegionName==region)]
        if len(pop)!=1: raise ValueError(f"need exactly one population for {country}/{region}")
        u = segment[list(included_ip)].ffill().fillna(0).to_numpy(float).T
        regression_start = 0 if regression_start_date is None else int(np.searchsorted(segment.Date.to_numpy(),np.datetime64(regression_start_date)))
        model = fit_npi_model(segment.ConfirmedCases.to_numpy(),u,float(pop.Population2020.iloc[0]),npi_maxes,regression_start=regression_start)
        model.update(country_name=country,region_name=region,last_date=segment.Date.iloc[-1].strftime('%Y-%m-%d'),last_input=u[:,-1].tolist())
        models.append(model)
    if not models: raise ValueError("no selected geography has at least 14 training days")
    bundle = dict(schema_version=1,npi_columns=list(included_ip),models=models)
    Path(trained_model_params_file).write_text(json.dumps(bundle,indent=2,allow_nan=False)+'\n')
    return bundle


def prescribe_npi(start_date_str,end_date_str,ip_file,costs_file,output_file, *, model_file,epsilon=.5):
    """Write XPRIZE-format daily prescriptions from a trained JSON bundle.

    ip_file supplies selected CountryName/RegionName pairs and historic
    interventions. costs_file has intervention weights in the bundle's
    column names. model_file is explicit, replacing historical hardcoded
    local paths. Returns the written DataFrame. No external service is used.
    """
    bundle = json.loads(Path(model_file).read_text()); costs = pd.read_csv(costs_file).fillna({'RegionName':''})
    plans = pd.read_csv(ip_file).fillna({'RegionName':''}); dates = pd.date_range(start_date_str,end_date_str)
    if dates.empty: raise ValueError("end date precedes start date")
    rows = []
    for model in bundle['models']:
        country,region = model['country_name'],model['region_name']
        if not ((plans.CountryName==country)&(plans.RegionName==region)).any(): continue
        if pd.Timestamp(model['last_date'])+pd.Timedelta(days=1)!=dates[0]:
            raise ValueError("prescription must begin the day after the training endpoint")
        selected = costs[(costs.CountryName==country)&(costs.RegionName==region)]
        if len(selected)!=1: raise ValueError(f"missing/duplicate costs for {country}/{region}")
        weights = selected[bundle['npi_columns']].iloc[0].to_numpy(float)
        control,_,_ = optimal_npi(model,len(dates),weights,epsilon)
        for k,date in enumerate(dates):
            rows.append(dict(CountryName=country,RegionName=region,Date=date.strftime('%Y-%m-%d'),**dict(zip(bundle['npi_columns'],control[:,k]))))
    if not rows: raise ValueError("no models matched the intervention plan")
    frame = pd.DataFrame(rows); frame.to_csv(output_file,index=False); return frame


def train_predict_prescribe_npi(npi_weights,human_npi_cost_factor,start_train_date_str,end_train_date_str,
        start_regression_date_str,end_predict_prescribe_date_str,data_file,geo_file,populations_file,
        included_ip,npi_mins,npi_maxes,trained_model_params_file):
    """Train and return per-geography fixed-policy and optimal-policy forecasts.

    Dates and CSV arguments follow train_npi_prescriptor. npi_mins/maxes are
    explicit bounds; output is a list of dicts with controls, cases, states,
    and (human,NPI) costs. All numerical work is in src, outside notebooks.
    """
    bundle = train_npi_prescriptor(start_train_date_str,end_train_date_str,data_file,geo_file,populations_file,
        included_ip,npi_maxes,trained_model_params_file,regression_start_date=start_regression_date_str)
    horizon = (pd.Timestamp(end_predict_prescribe_date_str)-pd.Timestamp(end_train_date_str)).days
    if horizon<1: raise ValueError("prediction endpoint must follow training endpoint")
    results = []
    for model in bundle['models']:
        model['params']['u_min'] = np.asarray(npi_mins,float).tolist()
        fixed = np.repeat(np.array(model['last_input'])[:,None],horizon,axis=1)
        cases,states = forecast_npi(model,fixed)
        control,opt_cases,opt_states = optimal_npi(model,horizon,npi_weights,human_npi_cost_factor)
        results.append(dict(country_name=model['country_name'],region_name=model['region_name'],fixed_inputs=fixed,
            fixed_cases=cases,fixed_states=states,optimal_inputs=control,optimal_cases=opt_cases,optimal_states=opt_states,
            fixed_cost=npi_cost(cases,fixed,npi_weights),optimal_cost=npi_cost(opt_cases,control,npi_weights)))
    return results


def forecast_quality_assessment(npi_weights,human_npi_cost_factor,start_train_date_str,end_train_date_str,
        start_regression_date_str,end_predict_prescribe_date_str,max_look_ahead_days,data_file,geo_file,
        populations_file,included_ip,npi_mins,npi_maxes,trained_model_params_file):
    """Evaluate held-out fixed-policy forecasts without fitting future observations.

    Returns a table with country/region, number of held-out days, MAE, RMSE,
    and bias. Actual held-out NPI schedules are known scenario inputs; case
    data after the training endpoint are used only for scoring.
    """
    bundle = train_npi_prescriptor(start_train_date_str,end_train_date_str,data_file,geo_file,populations_file,
        included_ip,npi_maxes,trained_model_params_file,regression_start_date=start_regression_date_str)
    data = read_oxford_data(data_file); rows = []
    end = min(pd.Timestamp(end_predict_prescribe_date_str),pd.Timestamp(end_train_date_str)+pd.Timedelta(days=max_look_ahead_days))
    for model in bundle['models']:
        segment = data[(data.CountryName==model['country_name'])&(data.RegionName==model['region_name'])&
            (data.Date>pd.Timestamp(end_train_date_str))&(data.Date<=end)]
        if segment.empty: continue
        if segment.Date.iloc[0]!=pd.Timestamp(end_train_date_str)+pd.Timedelta(days=1) or not (segment.Date.diff().dropna()==pd.Timedelta(days=1)).all():
            raise ValueError("held-out dates must be consecutive from the training endpoint")
        last = data[(data.CountryName==model['country_name'])&(data.RegionName==model['region_name'])&(data.Date==pd.Timestamp(end_train_date_str))].ConfirmedCases.iloc[0]
        actual = np.maximum(0,np.diff(np.r_[last,segment.ConfirmedCases.to_numpy()]))
        u = segment[list(included_ip)].ffill().fillna(0).to_numpy(float).T
        predicted,_ = forecast_npi(model,u); error = predicted-actual
        rows.append(dict(country_name=model['country_name'],region_name=model['region_name'],days=len(error),
                         mae=float(np.mean(abs(error))),rmse=float(np.sqrt(np.mean(error**2))),bias=float(np.mean(error))))
    return pd.DataFrame(rows)
