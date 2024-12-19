function res=Inv_params(params,monitor,varargin)
    if params.eq=="MT"
        params.ic=0.3;  
        params.num_per=20;
        params.to_cut=0.8;
        params.num_steps=500;
        params.reltol=6;
    else
        params.ic=0.01;  
        params.num_per=10;
        params.to_cut=0.8;
        params.num_steps=1500;
        params.reltol=6;
    end
    mdl=Train_model(params,monitor,varargin{:});
    res=mdl;
end

