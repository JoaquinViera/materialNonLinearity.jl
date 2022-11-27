var documenterSearchIndex = {"docs":
[{"location":"#FEA-with-general-material-models-on-beam-elements","page":"Home","title":"FEA with general material models on beam elements","text":"","category":"section"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Repository to perform finite element analyses of 2d noded Euler Bernoulli beam elements with linear and non-linear constitutive models.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [materialNonLinearity]","category":"page"},{"location":"#materialNonLinearity.AbstractAlgorithm","page":"Home","title":"materialNonLinearity.AbstractAlgorithm","text":"AbstractAlgorithm\n\nAbstract type to define numerical method algorithms.\n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.BoundaryConds","page":"Home","title":"materialNonLinearity.BoundaryConds","text":"BoundaryConds\n\nContains the boundary conditions of the structure.\nIncludes the support definition and the external forces.\n Fields:\n\nsuppMatrix::Matrix: Matrix whose i-th row contains the fixed degrees of freedom of the fixed i-th node according to the structure: [Nodeᵢ uz θy].  \nnodalForceMatrix::Matrix: Matrix whose i-th row contains the external nodal forces applied to the loaded i-th node according the structure: [Nodeᵢ Fz My].\n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.IsotropicBiLinear","page":"Home","title":"materialNonLinearity.IsotropicBiLinear","text":"IsotropicBiLinear\n\nIsotropic BiLinear constitutive model i.e. Linear softening or Linear Hardening. \nStress is proportional to strain given by either the Elasticity Modulus E and the Elastoplastic tangent modulus C.\nChange in proportionality occurs after the stress reaches a given value σY0.\n\nFields:\n\nE::Float64: Young Elasticity modulus.  \nσY0::Float64: Initial yield stress.\nK::Float64: Slope of the hardening function.\nne::Int64: Gauss integration points along the element length.\nns::Int64: Gauss integration points along the section height. \n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.IterParams","page":"Home","title":"materialNonLinearity.IterParams","text":"IterPArams\n\nStore the criterion stop at each load increment.\n\n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.LinearElastic","page":"Home","title":"materialNonLinearity.LinearElastic","text":"LinearElastic\n\nIsotropic Linear elastic constitutive model. \nStress is proportional to strain given by Elasticity Modulus E.\n\nFields:\n\nE::Float64: Young Elasticity modulus.  \nne::Int64: Gauss integration points along the element length.\nns::Int64: Gauss integration points along the section height. \n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.MaterialModel","page":"Home","title":"materialNonLinearity.MaterialModel","text":"MaterialModel\n\nAbstract type to define material constitutive model.\n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.Mesh","page":"Home","title":"materialNonLinearity.Mesh","text":"Mesh\n\nContains the mesh node coordinates and conectivity matrix.\n\n\nFields:\n\nnodesMat::Matrix{Float64}: Matrix containing the (xᵢ,zᵢ) coordinates of the i-th node. \nconecMat::Matrix: Matrix whose i-th row contains the information of the i-th element according the structure: [Material, Section, (nodeᵢ, nodeⱼ)].\n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.ModelSol","page":"Home","title":"materialNonLinearity.ModelSol","text":"ModelSol\n\nStores the solution of the problem.\n\n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.Rectangle","page":"Home","title":"materialNonLinearity.Rectangle","text":"Rectangle\n\nRectangular section.\n\n\nFields:\n\nb::Float64: Section width.\nh::Float64: Section height.\nA::Float64: Section area.\nIy::Float64: Second moment of inertia. \n\nNote:\nA and Iy are computed internally, only provide fields b and h.\n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.Section","page":"Home","title":"materialNonLinearity.Section","text":"Section\n\nAbstract type to define element sections.\n\n\n\n\n\n","category":"type"},{"location":"#materialNonLinearity.UserModel","page":"Home","title":"materialNonLinearity.UserModel","text":"UserModel\n\nUser defined constitutive model. \nThe user shall use import materialNonLinearity: constitutive_model and define the function constitutive_model as follows.\n\n\nfunction constitutive_model(ElemMaterialModel::UserModel, εₖ)\n    ...\n    return σ, ∂σ∂ε\nend\n\nFields:\n\nne::Int64: Gauss integration points along the element length.\nns::Int64: Gauss integration points along the section height. \n\n\n\n\n\n","category":"type"}]
}