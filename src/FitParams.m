classdef FitParams
    properties
        n;
        free;
        fix;
    end
    
    methods
        function obj = FitParams(n,varargin)
            %FITPARAMS Construct an instance of this class
            %   Detailed explanation goes here
            p=inputParser;
            p.addOptional("free",1:7);
            p.parse(varargin{:});
            free=p.Results.free;
            obj.n=n;
            obj.free=free;
            [~,pows]=Model.calculateMonomials(1,1,n);
            obj.fix=obj.InverseIndices(free,1:7);
            obj.free=[obj.InverseIndices(obj.fix,1:length(pows))];
        end
        function ind=InverseIndices(~,indices,full_indices)
            ind=full_indices;
            for i=full_indices
                if ismember(i,indices)
                    ind=ind(ind~=i);
                end
            end
        end
    end
end

