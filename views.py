from django.shortcuts import render, HttpResponse
from datetime import datetime
from home.models import Contact
from django.contrib import messages
from math import exp, sqrt
import numpy as np
from scipy.stats import norm
def index(request):
    if request.method == 'POST':
        S = float(request.POST['stock_price'])
        K = float(request.POST['strike_price'])
        r = float(request.POST['risk_free_rate'])
        sigma = float(request.POST['volatility'])
        T = float(request.POST['time'])
        N = int(request.POST['step_tree'])

        option_price = binomial_tree(S, K, r, sigma, T, N)

        return render(request, 'index.html', {'result': option_price})
    else:
        return render(request, 'index.html')

def n_step(request):
    if request.method == 'POST':
        S = float(request.POST['stock_price'])
        K = float(request.POST['strike_price'])
        r = float(request.POST['risk_free_rate'])
        sigma = float(request.POST['volatility'])
        T = float(request.POST['time'])
        N = int(request.POST['step_tree'])

        option_price = binomial_tree(S, K, r, sigma, T, N)

        return render(request, 'n_step.html', {'result': option_price})
    else:
        return render(request, 'n_step.html')
def about(request):
    return render(request, 'about.html') 
def base(request):
    return render (request, 'base.html')
def onestep(request):
    if request.method == 'POST':
        S = float(request.POST['stock_price'])
        K = float(request.POST['strike_price'])
        r = float(request.POST['risk_free_rate'])
        sigma = float(request.POST['volatility'])
        T = float(request.POST['time'])
        N = 1

        option_price = binomial_tree(S, K, r, sigma, T, N)

        return render(request, 'onestep.html', {'result': option_price})
    else:
        return render(request, 'onestep.html')

def twostep(request):
    if request.method == 'POST':
        S = float(request.POST['stock_price'])
        K = float(request.POST['strike_price'])
        r = float(request.POST['risk_free_rate'])
        sigma = float(request.POST['volatility'])
        T = float(request.POST['time'])
        N = 2

        option_price = binomial_tree(S, K, r, sigma, T, N)

        return render(request, 'twostep.html', {'result': option_price})
    else:
        return render(request, 'twostep.html')
    
def blackshole_call(request):
    if request.method == 'POST':
        S = float(request.POST['underlying_spot_price'])
        K = float(request.POST['strike_price'])
        T = float(request.POST['days_to_maturity'])
        sigma = float(request.POST['volatility'])
        r = float(request.POST['risk_free_rate'])
        option_price = blackshole_calculator_call(S, K, r, sigma, T)

        return render(request, 'blackshole_call.html', {'result':  option_price})
    else:
        return render(request, 'blackshole_call.html')

def blackshole_put(request):
    if request.method == 'POST':
        S = float(request.POST['underlying_spot_price'])
        K = float(request.POST['strike_price'])
        T = float(request.POST['days_to_maturity'])
        sigma = float(request.POST['volatility'])
        r = float(request.POST['risk_free_rate'])
        option_price = blackshole_calculator_put(S, K, r, sigma, T) 
        return render(request, 'blackshole_put.html', {'result': option_price})
    else:
        return render(request, 'blackshole_put.html',)
    
def contact(request):
    if request.method == "POST":
        name = request.POST.get('name')
        email = request.POST.get('email')
        phone = request.POST.get('phone')
        desc = request.POST.get('desc')
        contact = Contact(name=name, email=email, phone=phone, desc=desc, date = datetime.today())
        contact.save()
        messages.success(request, 'Your message has been sent!')
    return render(request, 'contact.html')
def blackshole_calculator_call(S,K,r,sigma,T):
     # cumulative function of standard normal distribution (risk-adjusted probability that the option will be exercised)     
     d1 = (np.log(S /K) + (r + 0.5 * sigma ** 2) *T) / (sigma * np.sqrt(T))

     # cumulative function of standard normal distribution (probability of receiving the stock at expiration of the option)
     d2 = (np.log(S /K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
     return (S * norm.cdf(d1, 0.0, 1.0) - K * np.exp(-r *T) * norm.cdf(d2, 0.0, 1.0))

def blackshole_calculator_put(S,K,r,sigma,T):
       # cumulative function of standard normal distribution (risk-adjusted probability that the option will be exercised)    
        d1 = (np.log(S /K) + (r + 0.5 *sigma ** 2) *T) / (sigma * np.sqrt(T))

        # cumulative function of standard normal distribution (probability of receiving the stock at expiration of the option)
        d2 = (np.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))

        return (K * np.exp(-r * T) * norm.cdf(-d2, 0.0, 1.0) - S * norm.cdf(-d1, 0.0, 1.0))

def binomial_tree(S, K, r, sigma, T, N):
    delta_t = T / N
    u = exp(sigma * sqrt(delta_t))
    d = 1 / u
    p = (exp(r * delta_t) - d) / (u - d)

    # Calculate stock prices at each node
    stock_prices = []
    for i in range(N + 1):
        stock_price = S * (u ** (N - i)) * (d ** i)
        stock_prices.append(stock_price)

    # Calculate option values at each node at time T
    option_values = []
    for i in range(N + 1):
        option_value = max(stock_prices[i] - K, 0)
        option_values.append(option_value)

    # Calculate option values at each node at time 0
    for n in range(N):
        for i in range(N - n):
            option_values[i] = exp(-r * delta_t) * (p * option_values[i] + (1 - p) * option_values[i + 1])

    return option_values[0]

