'''
Function to Calculate AIC

Bridgette Befort, Alexander Dowling, Ed Maginn (2021)
University of Notre Dame

'''

def AIC(SRS,n):
    
    '''
    Function to calculate AIC value using the equation:
        
        AIC = 2*n - 2*L, 
    
    where L is the maximum log likelihood value for the model. To calculate L:
    
    L = -n/2*ln(2*pi) - n/2*ln(sigma^2) - 1/(2*sigma^2)*SRS
   
    
    Input:
    SRS - sum of residuals squared (single value)
    n - number of parameters (single value)
    
    Output:
    AIC - Akaike Information Criteria Value
    '''
    
    #Calculate sigma^2
    sigma_sq = SRS/n
    
    #Calculate L
    L = -n/2*np.log(2*np.pi) - n/2*np.log(sigma_sq) - 1/(2*sigma_sq)*SRS
    
    # Calculate AIC
    AIC = 2*n - 2*L
    
    return AIC