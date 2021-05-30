function T=createOneHotEncoding(T,tableVariable)
    %
    % Code written by Christopher L. Stokely, January 30, 2019
    % Written in MATLAB R2018B.
    %
    % Command:
    % outputTable = createOneHotEncoding(T,tableVariable)
    % 
    % Input variable T needs to be a table and the tableVariable should be
    % a variable in that table.  tableVariable should be a variable that is
    % categorical but it does not have to be.  The code below converts the
    % variable to categorical if it is not already so.  A table will be 
    % returned that is the original input table without tableVariable, but
    % with new variables representing the one-hot encoded tableVariable.
    %
    % By one hot encoding, predictor importances can become very useful
    % when employing machine learning - from a model interpretability stand
    % -point. Being able to assign an importance to an individual category
    % can be useful and important in some cases.
    % 
    % For educational purposes, try looking into these Machine Learning
    % toolbox commands after building a model:
    % 1) oobPermutedPredictorImportance
    % 2) predictorImportance  (Be careful - this one is known to mislead)
    % 3) FeatureSelectionNCARRegression
    % 4) fsrnca or fscnca
    % 5) sequentialfs
    % 6) plotPartialDependence
    % 7) Individual Conditional Expectation (ICE) plots
    %
    % Note a MATLAB bug or oversight from MathWorks regarding having an 
    % underscore in the variable names that are in the table...
    % Note that the output table has new variables with labels that have an
    % underscore.  Removing these variables with "removevars" requires the
    % user to specify the column to be removed with the column number, not
    % the variable name.  Otherwise unintended columns will be deleted.
    %
    
    %%
    % determine if it is a table and throw an error if not
    if ~istable(T)
        error('Input table variable is not a table!')
    end
        
    %%
    % determine if the table variable to be encoded actually exists - throw an error otherwise
    if ~iscolumn(T.(genvarname(tableVariable)))
        error('Column variable does not exist in table!')
    end
    
    %%
    % do we want to make this a categorical variable?  how many instances are
    % there for this variable?
    numCategories=numel(unique(T.(genvarname(tableVariable))));

    fprintf('There are %i unique values in the tableVariable. \n', numCategories);

    % start creating one hot encoding variables if the number of categories after
    % variable expansion is reasonable
    
    if numCategories<100 % maximum of 99 categories for the variable splitting to be created
        
        tempCatVariable=T.(genvarname(tableVariable));
              
        % convert it to categorical in case it is not already
        tempCatVariable=categorical(tempCatVariable);
        
        %now get the unique categories that will be one hot encoded
        uniqueCategories=unique(tempCatVariable);

        %remove variable from table that is being one-hot encoded
        T=removevars(T,tableVariable);        
        
        % BEWARE - dynamic variables created using the EVAL command
        % Other useful commands to try out are genvarname,
        % matlab.lang.makeUniqueStrings, and matlab.lang.makeValidName
        for indexCategories=1:numel(uniqueCategories)
            oneHot=double((tempCatVariable==uniqueCategories(indexCategories))); %want a numeric value - not a Boolean value
            eval((strcat(char(tableVariable),'_',char(uniqueCategories(indexCategories)),'=oneHot;')));
            eval(strcat('T=addvars(T,',strcat(char(tableVariable),'_',char(uniqueCategories(indexCategories)),');')));            
        end
        % NOTE: the user cannot remove variables created that have an
        % underscore in their variable name
        
    else
        error('There are over 99 categories for this table variable.  Unable to proceed unless user increases conditional on line 42 to a more appropriate level.');
    end
        
end % end of function
