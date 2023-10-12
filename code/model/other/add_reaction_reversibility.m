% add revesibility into the model
% note: this kind of model is from deep learning and bayesian optimization


emodel = load('emodel_Saccharomyces_cerevisiae_Posterior_mean.mat')
model1 = emodel
model1.rev = zeros(length(model1.lb),1)
for i = 1: length(model1.lb)
    if model1.lb(i) == 0
        model1.rev(i) = 0
    else
        model1.rev(i) =1
    end
end

# save model


% other curation
% remove field
model1 = rmfield(model1, 'metNotes')
model1 = rmfield(model1, 'metMetaNetXID')



model1.metCharges  = transpose(model1.metCharges)
model1.metMetaNetXID  = transpose(model1.metMetaNetXID)
model1.rxnECNumbers  = transpose(model1.rxnECNumbers)
model1.rxnKEGGID  = transpose(model1.rxnKEGGID)
model1.rxnMetaNetXID  = transpose(model1.rxnMetaNetXID)
model1.excarbon  = transpose(model1.excarbon)


writeCbModel(model1, 'format','sbml', 'fileName', 'ecYeast.xml')


% for this kind of models, for gene infomation, do the following changes:
% replace __32__ as null
% replace YGL002C [c] as YGL002C



% error info
% metCharges: Y Size of metCharges was 3893. Expected 1
% metNotes: Y Size of metNotes was 3893. Expected 1
% metMetaNetXID: mets: Size of metMetaNetXID does not match elements in mets
% rxnECNumbers: Y Size of rxnECNumbers was 7359. Expected 1
% rxnKEGGID: Y Size of rxnKEGGID was 7359. Expected 1
% rxnMetaNetXID: Y Size of rxnMetaNetXID was 7359. Expected 1
