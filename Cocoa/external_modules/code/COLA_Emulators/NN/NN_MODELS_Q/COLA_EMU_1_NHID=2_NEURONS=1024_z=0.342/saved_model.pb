ݡ
��
^
AssignVariableOp
resource
value"dtype"
dtypetype"
validate_shapebool( �
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
N
Cast	
x"SrcT	
y"DstT"
SrcTtype"
DstTtype"
Truncatebool( 
8
Const
output"dtype"
valuetensor"
dtypetype
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(�

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
E
Relu
features"T
activations"T"
Ttype:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ��
@
StaticRegexFullMatch	
input

output
"
patternstring
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.8.02v2.8.0-rc1-32-g3f878cff5b68��

hid_layer1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*"
shared_namehid_layer1/kernel
x
%hid_layer1/kernel/Read/ReadVariableOpReadVariableOphid_layer1/kernel*
_output_shapes
:	�*
dtype0
w
hid_layer1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�* 
shared_namehid_layer1/bias
p
#hid_layer1/bias/Read/ReadVariableOpReadVariableOphid_layer1/bias*
_output_shapes	
:�*
dtype0
�
hid_layer2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*"
shared_namehid_layer2/kernel
y
%hid_layer2/kernel/Read/ReadVariableOpReadVariableOphid_layer2/kernel* 
_output_shapes
:
��*
dtype0
w
hid_layer2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:�* 
shared_namehid_layer2/bias
p
#hid_layer2/bias/Read/ReadVariableOpReadVariableOphid_layer2/bias*
_output_shapes	
:�*
dtype0
}
out_layer/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*!
shared_nameout_layer/kernel
v
$out_layer/kernel/Read/ReadVariableOpReadVariableOpout_layer/kernel*
_output_shapes
:	�*
dtype0
t
out_layer/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameout_layer/bias
m
"out_layer/bias/Read/ReadVariableOpReadVariableOpout_layer/bias*
_output_shapes
:*
dtype0
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
x
Adam/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/learning_rate
q
&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*
_output_shapes
: *
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
�
Adam/hid_layer1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*)
shared_nameAdam/hid_layer1/kernel/m
�
,Adam/hid_layer1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/hid_layer1/kernel/m*
_output_shapes
:	�*
dtype0
�
Adam/hid_layer1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*'
shared_nameAdam/hid_layer1/bias/m
~
*Adam/hid_layer1/bias/m/Read/ReadVariableOpReadVariableOpAdam/hid_layer1/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/hid_layer2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*)
shared_nameAdam/hid_layer2/kernel/m
�
,Adam/hid_layer2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/hid_layer2/kernel/m* 
_output_shapes
:
��*
dtype0
�
Adam/hid_layer2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*'
shared_nameAdam/hid_layer2/bias/m
~
*Adam/hid_layer2/bias/m/Read/ReadVariableOpReadVariableOpAdam/hid_layer2/bias/m*
_output_shapes	
:�*
dtype0
�
Adam/out_layer/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*(
shared_nameAdam/out_layer/kernel/m
�
+Adam/out_layer/kernel/m/Read/ReadVariableOpReadVariableOpAdam/out_layer/kernel/m*
_output_shapes
:	�*
dtype0
�
Adam/out_layer/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*&
shared_nameAdam/out_layer/bias/m
{
)Adam/out_layer/bias/m/Read/ReadVariableOpReadVariableOpAdam/out_layer/bias/m*
_output_shapes
:*
dtype0
�
Adam/hid_layer1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*)
shared_nameAdam/hid_layer1/kernel/v
�
,Adam/hid_layer1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/hid_layer1/kernel/v*
_output_shapes
:	�*
dtype0
�
Adam/hid_layer1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*'
shared_nameAdam/hid_layer1/bias/v
~
*Adam/hid_layer1/bias/v/Read/ReadVariableOpReadVariableOpAdam/hid_layer1/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/hid_layer2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
��*)
shared_nameAdam/hid_layer2/kernel/v
�
,Adam/hid_layer2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/hid_layer2/kernel/v* 
_output_shapes
:
��*
dtype0
�
Adam/hid_layer2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:�*'
shared_nameAdam/hid_layer2/bias/v
~
*Adam/hid_layer2/bias/v/Read/ReadVariableOpReadVariableOpAdam/hid_layer2/bias/v*
_output_shapes	
:�*
dtype0
�
Adam/out_layer/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*(
shared_nameAdam/out_layer/kernel/v
�
+Adam/out_layer/kernel/v/Read/ReadVariableOpReadVariableOpAdam/out_layer/kernel/v*
_output_shapes
:	�*
dtype0
�
Adam/out_layer/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*&
shared_nameAdam/out_layer/bias/v
{
)Adam/out_layer/bias/v/Read/ReadVariableOpReadVariableOpAdam/out_layer/bias/v*
_output_shapes
:*
dtype0

NoOpNoOp
�+
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�*
value�*B�* B�*
�
layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
	optimizer
	variables
trainable_variables
regularization_losses
	keras_api
	__call__
*
&call_and_return_all_conditional_losses
_default_save_signature

signatures*
�

kernel
bias
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses*
�

kernel
bias
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses*
�

kernel
bias
	variables
 trainable_variables
!regularization_losses
"	keras_api
#__call__
*$&call_and_return_all_conditional_losses*
�
%iter

&beta_1

'beta_2
	(decay
)learning_ratemDmEmFmGmHmIvJvKvLvMvNvO*
.
0
1
2
3
4
5*
.
0
1
2
3
4
5*
* 
�
*non_trainable_variables

+layers
,metrics
-layer_regularization_losses
.layer_metrics
	variables
trainable_variables
regularization_losses
	__call__
_default_save_signature
*
&call_and_return_all_conditional_losses
&
"call_and_return_conditional_losses*
* 
* 
* 

/serving_default* 
a[
VARIABLE_VALUEhid_layer1/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEhid_layer1/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE*

0
1*

0
1*
* 
�
0non_trainable_variables

1layers
2metrics
3layer_regularization_losses
4layer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*
* 
* 
a[
VARIABLE_VALUEhid_layer2/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEhid_layer2/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE*

0
1*

0
1*
* 
�
5non_trainable_variables

6layers
7metrics
8layer_regularization_losses
9layer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*
* 
* 
`Z
VARIABLE_VALUEout_layer/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEout_layer/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE*

0
1*

0
1*
* 
�
:non_trainable_variables

;layers
<metrics
=layer_regularization_losses
>layer_metrics
	variables
 trainable_variables
!regularization_losses
#__call__
*$&call_and_return_all_conditional_losses
&$"call_and_return_conditional_losses*
* 
* 
LF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE*
PJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE*
PJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE*
NH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
* 

0
1
2*

?0*
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
8
	@total
	Acount
B	variables
C	keras_api*
SM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
SM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

@0
A1*

B	variables*
�~
VARIABLE_VALUEAdam/hid_layer1/kernel/mRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUEAdam/hid_layer1/bias/mPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
�~
VARIABLE_VALUEAdam/hid_layer2/kernel/mRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUEAdam/hid_layer2/bias/mPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
�}
VARIABLE_VALUEAdam/out_layer/kernel/mRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUEAdam/out_layer/bias/mPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
�~
VARIABLE_VALUEAdam/hid_layer1/kernel/vRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUEAdam/hid_layer1/bias/vPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
�~
VARIABLE_VALUEAdam/hid_layer2/kernel/vRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
�z
VARIABLE_VALUEAdam/hid_layer2/bias/vPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
�}
VARIABLE_VALUEAdam/out_layer/kernel/vRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
y
VARIABLE_VALUEAdam/out_layer/bias/vPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
�
 serving_default_hid_layer1_inputPlaceholder*'
_output_shapes
:���������*
dtype0*
shape:���������
�
StatefulPartitionedCallStatefulPartitionedCall serving_default_hid_layer1_inputhid_layer1/kernelhid_layer1/biashid_layer2/kernelhid_layer2/biasout_layer/kernelout_layer/bias*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *.
f)R'
%__inference_signature_wrapper_3340493
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�

StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename%hid_layer1/kernel/Read/ReadVariableOp#hid_layer1/bias/Read/ReadVariableOp%hid_layer2/kernel/Read/ReadVariableOp#hid_layer2/bias/Read/ReadVariableOp$out_layer/kernel/Read/ReadVariableOp"out_layer/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOp,Adam/hid_layer1/kernel/m/Read/ReadVariableOp*Adam/hid_layer1/bias/m/Read/ReadVariableOp,Adam/hid_layer2/kernel/m/Read/ReadVariableOp*Adam/hid_layer2/bias/m/Read/ReadVariableOp+Adam/out_layer/kernel/m/Read/ReadVariableOp)Adam/out_layer/bias/m/Read/ReadVariableOp,Adam/hid_layer1/kernel/v/Read/ReadVariableOp*Adam/hid_layer1/bias/v/Read/ReadVariableOp,Adam/hid_layer2/kernel/v/Read/ReadVariableOp*Adam/hid_layer2/bias/v/Read/ReadVariableOp+Adam/out_layer/kernel/v/Read/ReadVariableOp)Adam/out_layer/bias/v/Read/ReadVariableOpConst*&
Tin
2	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *)
f$R"
 __inference__traced_save_3340651
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamehid_layer1/kernelhid_layer1/biashid_layer2/kernelhid_layer2/biasout_layer/kernelout_layer/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_ratetotalcountAdam/hid_layer1/kernel/mAdam/hid_layer1/bias/mAdam/hid_layer2/kernel/mAdam/hid_layer2/bias/mAdam/out_layer/kernel/mAdam/out_layer/bias/mAdam/hid_layer1/kernel/vAdam/hid_layer1/bias/vAdam/hid_layer2/kernel/vAdam/hid_layer2/bias/vAdam/out_layer/kernel/vAdam/out_layer/bias/v*%
Tin
2*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *,
f'R%
#__inference__traced_restore_3340736��
�%
�
"__inference__wrapped_model_3340166
hid_layer1_inputJ
7sequential_17_hid_layer1_matmul_readvariableop_resource:	�G
8sequential_17_hid_layer1_biasadd_readvariableop_resource:	�K
7sequential_17_hid_layer2_matmul_readvariableop_resource:
��G
8sequential_17_hid_layer2_biasadd_readvariableop_resource:	�I
6sequential_17_out_layer_matmul_readvariableop_resource:	�E
7sequential_17_out_layer_biasadd_readvariableop_resource:
identity��/sequential_17/hid_layer1/BiasAdd/ReadVariableOp�.sequential_17/hid_layer1/MatMul/ReadVariableOp�/sequential_17/hid_layer2/BiasAdd/ReadVariableOp�.sequential_17/hid_layer2/MatMul/ReadVariableOp�.sequential_17/out_layer/BiasAdd/ReadVariableOp�-sequential_17/out_layer/MatMul/ReadVariableOpm
sequential_17/CastCasthid_layer1_input*

DstT0*

SrcT0*'
_output_shapes
:����������
.sequential_17/hid_layer1/MatMul/ReadVariableOpReadVariableOp7sequential_17_hid_layer1_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_17/hid_layer1/MatMulMatMulsequential_17/Cast:y:06sequential_17/hid_layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
/sequential_17/hid_layer1/BiasAdd/ReadVariableOpReadVariableOp8sequential_17_hid_layer1_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 sequential_17/hid_layer1/BiasAddBiasAdd)sequential_17/hid_layer1/MatMul:product:07sequential_17/hid_layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_17/hid_layer1/ReluRelu)sequential_17/hid_layer1/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
.sequential_17/hid_layer2/MatMul/ReadVariableOpReadVariableOp7sequential_17_hid_layer2_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
sequential_17/hid_layer2/MatMulMatMul+sequential_17/hid_layer1/Relu:activations:06sequential_17/hid_layer2/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
/sequential_17/hid_layer2/BiasAdd/ReadVariableOpReadVariableOp8sequential_17_hid_layer2_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
 sequential_17/hid_layer2/BiasAddBiasAdd)sequential_17/hid_layer2/MatMul:product:07sequential_17/hid_layer2/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
sequential_17/hid_layer2/ReluRelu)sequential_17/hid_layer2/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
-sequential_17/out_layer/MatMul/ReadVariableOpReadVariableOp6sequential_17_out_layer_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
sequential_17/out_layer/MatMulMatMul+sequential_17/hid_layer2/Relu:activations:05sequential_17/out_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
.sequential_17/out_layer/BiasAdd/ReadVariableOpReadVariableOp7sequential_17_out_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
sequential_17/out_layer/BiasAddBiasAdd(sequential_17/out_layer/MatMul:product:06sequential_17/out_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
sequential_17/out_layer/ReluRelu(sequential_17/out_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������y
IdentityIdentity*sequential_17/out_layer/Relu:activations:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp0^sequential_17/hid_layer1/BiasAdd/ReadVariableOp/^sequential_17/hid_layer1/MatMul/ReadVariableOp0^sequential_17/hid_layer2/BiasAdd/ReadVariableOp/^sequential_17/hid_layer2/MatMul/ReadVariableOp/^sequential_17/out_layer/BiasAdd/ReadVariableOp.^sequential_17/out_layer/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 2b
/sequential_17/hid_layer1/BiasAdd/ReadVariableOp/sequential_17/hid_layer1/BiasAdd/ReadVariableOp2`
.sequential_17/hid_layer1/MatMul/ReadVariableOp.sequential_17/hid_layer1/MatMul/ReadVariableOp2b
/sequential_17/hid_layer2/BiasAdd/ReadVariableOp/sequential_17/hid_layer2/BiasAdd/ReadVariableOp2`
.sequential_17/hid_layer2/MatMul/ReadVariableOp.sequential_17/hid_layer2/MatMul/ReadVariableOp2`
.sequential_17/out_layer/BiasAdd/ReadVariableOp.sequential_17/out_layer/BiasAdd/ReadVariableOp2^
-sequential_17/out_layer/MatMul/ReadVariableOp-sequential_17/out_layer/MatMul/ReadVariableOp:Y U
'
_output_shapes
:���������
*
_user_specified_namehid_layer1_input
�
�
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340226

inputs%
hid_layer1_3340186:	�!
hid_layer1_3340188:	�&
hid_layer2_3340203:
��!
hid_layer2_3340205:	�$
out_layer_3340220:	�
out_layer_3340222:
identity��"hid_layer1/StatefulPartitionedCall�"hid_layer2/StatefulPartitionedCall�!out_layer/StatefulPartitionedCallU
CastCastinputs*

DstT0*

SrcT0*'
_output_shapes
:����������
"hid_layer1/StatefulPartitionedCallStatefulPartitionedCallCast:y:0hid_layer1_3340186hid_layer1_3340188*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer1_layer_call_and_return_conditional_losses_3340185�
"hid_layer2/StatefulPartitionedCallStatefulPartitionedCall+hid_layer1/StatefulPartitionedCall:output:0hid_layer2_3340203hid_layer2_3340205*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer2_layer_call_and_return_conditional_losses_3340202�
!out_layer/StatefulPartitionedCallStatefulPartitionedCall+hid_layer2/StatefulPartitionedCall:output:0out_layer_3340220out_layer_3340222*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_out_layer_layer_call_and_return_conditional_losses_3340219y
IdentityIdentity*out_layer/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp#^hid_layer1/StatefulPartitionedCall#^hid_layer2/StatefulPartitionedCall"^out_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 2H
"hid_layer1/StatefulPartitionedCall"hid_layer1/StatefulPartitionedCall2H
"hid_layer2/StatefulPartitionedCall"hid_layer2/StatefulPartitionedCall2F
!out_layer/StatefulPartitionedCall!out_layer/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340310

inputs%
hid_layer1_3340294:	�!
hid_layer1_3340296:	�&
hid_layer2_3340299:
��!
hid_layer2_3340301:	�$
out_layer_3340304:	�
out_layer_3340306:
identity��"hid_layer1/StatefulPartitionedCall�"hid_layer2/StatefulPartitionedCall�!out_layer/StatefulPartitionedCallU
CastCastinputs*

DstT0*

SrcT0*'
_output_shapes
:����������
"hid_layer1/StatefulPartitionedCallStatefulPartitionedCallCast:y:0hid_layer1_3340294hid_layer1_3340296*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer1_layer_call_and_return_conditional_losses_3340185�
"hid_layer2/StatefulPartitionedCallStatefulPartitionedCall+hid_layer1/StatefulPartitionedCall:output:0hid_layer2_3340299hid_layer2_3340301*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer2_layer_call_and_return_conditional_losses_3340202�
!out_layer/StatefulPartitionedCallStatefulPartitionedCall+hid_layer2/StatefulPartitionedCall:output:0out_layer_3340304out_layer_3340306*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_out_layer_layer_call_and_return_conditional_losses_3340219y
IdentityIdentity*out_layer/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp#^hid_layer1/StatefulPartitionedCall#^hid_layer2/StatefulPartitionedCall"^out_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 2H
"hid_layer1/StatefulPartitionedCall"hid_layer1/StatefulPartitionedCall2H
"hid_layer2/StatefulPartitionedCall"hid_layer2/StatefulPartitionedCall2F
!out_layer/StatefulPartitionedCall!out_layer/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
%__inference_signature_wrapper_3340493
hid_layer1_input
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:	�
	unknown_4:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallhid_layer1_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *+
f&R$
"__inference__wrapped_model_3340166o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Y U
'
_output_shapes
:���������
*
_user_specified_namehid_layer1_input
�
�
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340448

inputs<
)hid_layer1_matmul_readvariableop_resource:	�9
*hid_layer1_biasadd_readvariableop_resource:	�=
)hid_layer2_matmul_readvariableop_resource:
��9
*hid_layer2_biasadd_readvariableop_resource:	�;
(out_layer_matmul_readvariableop_resource:	�7
)out_layer_biasadd_readvariableop_resource:
identity��!hid_layer1/BiasAdd/ReadVariableOp� hid_layer1/MatMul/ReadVariableOp�!hid_layer2/BiasAdd/ReadVariableOp� hid_layer2/MatMul/ReadVariableOp� out_layer/BiasAdd/ReadVariableOp�out_layer/MatMul/ReadVariableOpU
CastCastinputs*

DstT0*

SrcT0*'
_output_shapes
:����������
 hid_layer1/MatMul/ReadVariableOpReadVariableOp)hid_layer1_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
hid_layer1/MatMulMatMulCast:y:0(hid_layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
!hid_layer1/BiasAdd/ReadVariableOpReadVariableOp*hid_layer1_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
hid_layer1/BiasAddBiasAddhid_layer1/MatMul:product:0)hid_layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������g
hid_layer1/ReluReluhid_layer1/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
 hid_layer2/MatMul/ReadVariableOpReadVariableOp)hid_layer2_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
hid_layer2/MatMulMatMulhid_layer1/Relu:activations:0(hid_layer2/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
!hid_layer2/BiasAdd/ReadVariableOpReadVariableOp*hid_layer2_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
hid_layer2/BiasAddBiasAddhid_layer2/MatMul:product:0)hid_layer2/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������g
hid_layer2/ReluReluhid_layer2/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
out_layer/MatMul/ReadVariableOpReadVariableOp(out_layer_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
out_layer/MatMulMatMulhid_layer2/Relu:activations:0'out_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 out_layer/BiasAdd/ReadVariableOpReadVariableOp)out_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
out_layer/BiasAddBiasAddout_layer/MatMul:product:0(out_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������d
out_layer/ReluReluout_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������k
IdentityIdentityout_layer/Relu:activations:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp"^hid_layer1/BiasAdd/ReadVariableOp!^hid_layer1/MatMul/ReadVariableOp"^hid_layer2/BiasAdd/ReadVariableOp!^hid_layer2/MatMul/ReadVariableOp!^out_layer/BiasAdd/ReadVariableOp ^out_layer/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 2F
!hid_layer1/BiasAdd/ReadVariableOp!hid_layer1/BiasAdd/ReadVariableOp2D
 hid_layer1/MatMul/ReadVariableOp hid_layer1/MatMul/ReadVariableOp2F
!hid_layer2/BiasAdd/ReadVariableOp!hid_layer2/BiasAdd/ReadVariableOp2D
 hid_layer2/MatMul/ReadVariableOp hid_layer2/MatMul/ReadVariableOp2D
 out_layer/BiasAdd/ReadVariableOp out_layer/BiasAdd/ReadVariableOp2B
out_layer/MatMul/ReadVariableOpout_layer/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
G__inference_hid_layer2_layer_call_and_return_conditional_losses_3340533

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:����������b
IdentityIdentityRelu:activations:0^NoOp*
T0*(
_output_shapes
:����������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�9
�

 __inference__traced_save_3340651
file_prefix0
,savev2_hid_layer1_kernel_read_readvariableop.
*savev2_hid_layer1_bias_read_readvariableop0
,savev2_hid_layer2_kernel_read_readvariableop.
*savev2_hid_layer2_bias_read_readvariableop/
+savev2_out_layer_kernel_read_readvariableop-
)savev2_out_layer_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop7
3savev2_adam_hid_layer1_kernel_m_read_readvariableop5
1savev2_adam_hid_layer1_bias_m_read_readvariableop7
3savev2_adam_hid_layer2_kernel_m_read_readvariableop5
1savev2_adam_hid_layer2_bias_m_read_readvariableop6
2savev2_adam_out_layer_kernel_m_read_readvariableop4
0savev2_adam_out_layer_bias_m_read_readvariableop7
3savev2_adam_hid_layer1_kernel_v_read_readvariableop5
1savev2_adam_hid_layer1_bias_v_read_readvariableop7
3savev2_adam_hid_layer2_kernel_v_read_readvariableop5
1savev2_adam_hid_layer2_bias_v_read_readvariableop6
2savev2_adam_out_layer_kernel_v_read_readvariableop4
0savev2_adam_out_layer_bias_v_read_readvariableop
savev2_const

identity_1��MergeV2Checkpointsw
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*Z
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.parta
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: f

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: L

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :f
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : �
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*�
value�B�B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*G
value>B<B B B B B B B B B B B B B B B B B B B B B B B B B B �

SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0,savev2_hid_layer1_kernel_read_readvariableop*savev2_hid_layer1_bias_read_readvariableop,savev2_hid_layer2_kernel_read_readvariableop*savev2_hid_layer2_bias_read_readvariableop+savev2_out_layer_kernel_read_readvariableop)savev2_out_layer_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop3savev2_adam_hid_layer1_kernel_m_read_readvariableop1savev2_adam_hid_layer1_bias_m_read_readvariableop3savev2_adam_hid_layer2_kernel_m_read_readvariableop1savev2_adam_hid_layer2_bias_m_read_readvariableop2savev2_adam_out_layer_kernel_m_read_readvariableop0savev2_adam_out_layer_bias_m_read_readvariableop3savev2_adam_hid_layer1_kernel_v_read_readvariableop1savev2_adam_hid_layer1_bias_v_read_readvariableop3savev2_adam_hid_layer2_kernel_v_read_readvariableop1savev2_adam_hid_layer2_bias_v_read_readvariableop2savev2_adam_out_layer_kernel_v_read_readvariableop0savev2_adam_out_layer_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *(
dtypes
2	�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 f
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: Q

Identity_1IdentityIdentity:output:0^NoOp*
T0*
_output_shapes
: [
NoOpNoOp^MergeV2Checkpoints*"
_acd_function_control_output(*
_output_shapes
 "!

identity_1Identity_1:output:0*�
_input_shapes�
�: :	�:�:
��:�:	�:: : : : : : : :	�:�:
��:�:	�::	�:�:
��:�:	�:: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:%!

_output_shapes
:	�:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:%!

_output_shapes
:	�: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :%!

_output_shapes
:	�:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:%!

_output_shapes
:	�: 

_output_shapes
::%!

_output_shapes
:	�:!

_output_shapes	
:�:&"
 
_output_shapes
:
��:!

_output_shapes	
:�:%!

_output_shapes
:	�: 

_output_shapes
::

_output_shapes
: 
�
�
,__inference_hid_layer1_layer_call_fn_3340502

inputs
unknown:	�
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer1_layer_call_and_return_conditional_losses_3340185p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
/__inference_sequential_17_layer_call_fn_3340342
hid_layer1_input
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:	�
	unknown_4:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallhid_layer1_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340310o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Y U
'
_output_shapes
:���������
*
_user_specified_namehid_layer1_input
�f
�
#__inference__traced_restore_3340736
file_prefix5
"assignvariableop_hid_layer1_kernel:	�1
"assignvariableop_1_hid_layer1_bias:	�8
$assignvariableop_2_hid_layer2_kernel:
��1
"assignvariableop_3_hid_layer2_bias:	�6
#assignvariableop_4_out_layer_kernel:	�/
!assignvariableop_5_out_layer_bias:&
assignvariableop_6_adam_iter:	 (
assignvariableop_7_adam_beta_1: (
assignvariableop_8_adam_beta_2: '
assignvariableop_9_adam_decay: 0
&assignvariableop_10_adam_learning_rate: #
assignvariableop_11_total: #
assignvariableop_12_count: ?
,assignvariableop_13_adam_hid_layer1_kernel_m:	�9
*assignvariableop_14_adam_hid_layer1_bias_m:	�@
,assignvariableop_15_adam_hid_layer2_kernel_m:
��9
*assignvariableop_16_adam_hid_layer2_bias_m:	�>
+assignvariableop_17_adam_out_layer_kernel_m:	�7
)assignvariableop_18_adam_out_layer_bias_m:?
,assignvariableop_19_adam_hid_layer1_kernel_v:	�9
*assignvariableop_20_adam_hid_layer1_bias_v:	�@
,assignvariableop_21_adam_hid_layer2_kernel_v:
��9
*assignvariableop_22_adam_hid_layer2_bias_v:	�>
+assignvariableop_23_adam_out_layer_kernel_v:	�7
)assignvariableop_24_adam_out_layer_bias_v:
identity_26��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_3�AssignVariableOp_4�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*�
value�B�B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*G
value>B<B B B B B B B B B B B B B B B B B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*|
_output_shapesj
h::::::::::::::::::::::::::*(
dtypes
2	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOp"assignvariableop_hid_layer1_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOp"assignvariableop_1_hid_layer1_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOp$assignvariableop_2_hid_layer2_kernelIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOp"assignvariableop_3_hid_layer2_biasIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOp#assignvariableop_4_out_layer_kernelIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOp!assignvariableop_5_out_layer_biasIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOp_6AssignVariableOpassignvariableop_6_adam_iterIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOpassignvariableop_7_adam_beta_1Identity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOpassignvariableop_8_adam_beta_2Identity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOpassignvariableop_9_adam_decayIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOp&assignvariableop_10_adam_learning_rateIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOpassignvariableop_11_totalIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOpassignvariableop_12_countIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOp,assignvariableop_13_adam_hid_layer1_kernel_mIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_14AssignVariableOp*assignvariableop_14_adam_hid_layer1_bias_mIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOp,assignvariableop_15_adam_hid_layer2_kernel_mIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOp*assignvariableop_16_adam_hid_layer2_bias_mIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOp+assignvariableop_17_adam_out_layer_kernel_mIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOp)assignvariableop_18_adam_out_layer_bias_mIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOp,assignvariableop_19_adam_hid_layer1_kernel_vIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOp*assignvariableop_20_adam_hid_layer1_bias_vIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOp,assignvariableop_21_adam_hid_layer2_kernel_vIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOp*assignvariableop_22_adam_hid_layer2_bias_vIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOp+assignvariableop_23_adam_out_layer_kernel_vIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOp)assignvariableop_24_adam_out_layer_bias_vIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype01
NoOpNoOp"/device:CPU:0*
_output_shapes
 �
Identity_25Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_26IdentityIdentity_25:output:0^NoOp_1*
T0*
_output_shapes
: �
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 "#
identity_26Identity_26:output:0*G
_input_shapes6
4: : : : : : : : : : : : : : : : : : : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�

�
F__inference_out_layer_layer_call_and_return_conditional_losses_3340553

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
+__inference_out_layer_layer_call_fn_3340542

inputs
unknown:	�
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_out_layer_layer_call_and_return_conditional_losses_3340219o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
/__inference_sequential_17_layer_call_fn_3340405

inputs
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:	�
	unknown_4:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340226o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
,__inference_hid_layer2_layer_call_fn_3340522

inputs
unknown:
��
	unknown_0:	�
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer2_layer_call_and_return_conditional_losses_3340202p
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:����������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340362
hid_layer1_input%
hid_layer1_3340346:	�!
hid_layer1_3340348:	�&
hid_layer2_3340351:
��!
hid_layer2_3340353:	�$
out_layer_3340356:	�
out_layer_3340358:
identity��"hid_layer1/StatefulPartitionedCall�"hid_layer2/StatefulPartitionedCall�!out_layer/StatefulPartitionedCall_
CastCasthid_layer1_input*

DstT0*

SrcT0*'
_output_shapes
:����������
"hid_layer1/StatefulPartitionedCallStatefulPartitionedCallCast:y:0hid_layer1_3340346hid_layer1_3340348*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer1_layer_call_and_return_conditional_losses_3340185�
"hid_layer2/StatefulPartitionedCallStatefulPartitionedCall+hid_layer1/StatefulPartitionedCall:output:0hid_layer2_3340351hid_layer2_3340353*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer2_layer_call_and_return_conditional_losses_3340202�
!out_layer/StatefulPartitionedCallStatefulPartitionedCall+hid_layer2/StatefulPartitionedCall:output:0out_layer_3340356out_layer_3340358*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_out_layer_layer_call_and_return_conditional_losses_3340219y
IdentityIdentity*out_layer/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp#^hid_layer1/StatefulPartitionedCall#^hid_layer2/StatefulPartitionedCall"^out_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 2H
"hid_layer1/StatefulPartitionedCall"hid_layer1/StatefulPartitionedCall2H
"hid_layer2/StatefulPartitionedCall"hid_layer2/StatefulPartitionedCall2F
!out_layer/StatefulPartitionedCall!out_layer/StatefulPartitionedCall:Y U
'
_output_shapes
:���������
*
_user_specified_namehid_layer1_input
�

�
G__inference_hid_layer1_layer_call_and_return_conditional_losses_3340185

inputs1
matmul_readvariableop_resource:	�.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:����������b
IdentityIdentityRelu:activations:0^NoOp*
T0*(
_output_shapes
:����������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
G__inference_hid_layer2_layer_call_and_return_conditional_losses_3340202

inputs2
matmul_readvariableop_resource:
��.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpv
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:����������b
IdentityIdentityRelu:activations:0^NoOp*
T0*(
_output_shapes
:����������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�

�
F__inference_out_layer_layer_call_and_return_conditional_losses_3340219

inputs1
matmul_readvariableop_resource:	�-
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:����������
 
_user_specified_nameinputs
�
�
/__inference_sequential_17_layer_call_fn_3340422

inputs
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:	�
	unknown_4:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340310o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�	
�
/__inference_sequential_17_layer_call_fn_3340241
hid_layer1_input
unknown:	�
	unknown_0:	�
	unknown_1:
��
	unknown_2:	�
	unknown_3:	�
	unknown_4:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallhid_layer1_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *S
fNRL
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340226o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:Y U
'
_output_shapes
:���������
*
_user_specified_namehid_layer1_input
�
�
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340474

inputs<
)hid_layer1_matmul_readvariableop_resource:	�9
*hid_layer1_biasadd_readvariableop_resource:	�=
)hid_layer2_matmul_readvariableop_resource:
��9
*hid_layer2_biasadd_readvariableop_resource:	�;
(out_layer_matmul_readvariableop_resource:	�7
)out_layer_biasadd_readvariableop_resource:
identity��!hid_layer1/BiasAdd/ReadVariableOp� hid_layer1/MatMul/ReadVariableOp�!hid_layer2/BiasAdd/ReadVariableOp� hid_layer2/MatMul/ReadVariableOp� out_layer/BiasAdd/ReadVariableOp�out_layer/MatMul/ReadVariableOpU
CastCastinputs*

DstT0*

SrcT0*'
_output_shapes
:����������
 hid_layer1/MatMul/ReadVariableOpReadVariableOp)hid_layer1_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
hid_layer1/MatMulMatMulCast:y:0(hid_layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
!hid_layer1/BiasAdd/ReadVariableOpReadVariableOp*hid_layer1_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
hid_layer1/BiasAddBiasAddhid_layer1/MatMul:product:0)hid_layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������g
hid_layer1/ReluReluhid_layer1/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
 hid_layer2/MatMul/ReadVariableOpReadVariableOp)hid_layer2_matmul_readvariableop_resource* 
_output_shapes
:
��*
dtype0�
hid_layer2/MatMulMatMulhid_layer1/Relu:activations:0(hid_layer2/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
!hid_layer2/BiasAdd/ReadVariableOpReadVariableOp*hid_layer2_biasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0�
hid_layer2/BiasAddBiasAddhid_layer2/MatMul:product:0)hid_layer2/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������g
hid_layer2/ReluReluhid_layer2/BiasAdd:output:0*
T0*(
_output_shapes
:�����������
out_layer/MatMul/ReadVariableOpReadVariableOp(out_layer_matmul_readvariableop_resource*
_output_shapes
:	�*
dtype0�
out_layer/MatMulMatMulhid_layer2/Relu:activations:0'out_layer/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 out_layer/BiasAdd/ReadVariableOpReadVariableOp)out_layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
out_layer/BiasAddBiasAddout_layer/MatMul:product:0(out_layer/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������d
out_layer/ReluReluout_layer/BiasAdd:output:0*
T0*'
_output_shapes
:���������k
IdentityIdentityout_layer/Relu:activations:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp"^hid_layer1/BiasAdd/ReadVariableOp!^hid_layer1/MatMul/ReadVariableOp"^hid_layer2/BiasAdd/ReadVariableOp!^hid_layer2/MatMul/ReadVariableOp!^out_layer/BiasAdd/ReadVariableOp ^out_layer/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 2F
!hid_layer1/BiasAdd/ReadVariableOp!hid_layer1/BiasAdd/ReadVariableOp2D
 hid_layer1/MatMul/ReadVariableOp hid_layer1/MatMul/ReadVariableOp2F
!hid_layer2/BiasAdd/ReadVariableOp!hid_layer2/BiasAdd/ReadVariableOp2D
 hid_layer2/MatMul/ReadVariableOp hid_layer2/MatMul/ReadVariableOp2D
 out_layer/BiasAdd/ReadVariableOp out_layer/BiasAdd/ReadVariableOp2B
out_layer/MatMul/ReadVariableOpout_layer/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
G__inference_hid_layer1_layer_call_and_return_conditional_losses_3340513

inputs1
matmul_readvariableop_resource:	�.
biasadd_readvariableop_resource:	�
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�*
dtype0j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������s
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:�*
dtype0w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������Q
ReluReluBiasAdd:output:0*
T0*(
_output_shapes
:����������b
IdentityIdentityRelu:activations:0^NoOp*
T0*(
_output_shapes
:����������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340382
hid_layer1_input%
hid_layer1_3340366:	�!
hid_layer1_3340368:	�&
hid_layer2_3340371:
��!
hid_layer2_3340373:	�$
out_layer_3340376:	�
out_layer_3340378:
identity��"hid_layer1/StatefulPartitionedCall�"hid_layer2/StatefulPartitionedCall�!out_layer/StatefulPartitionedCall_
CastCasthid_layer1_input*

DstT0*

SrcT0*'
_output_shapes
:����������
"hid_layer1/StatefulPartitionedCallStatefulPartitionedCallCast:y:0hid_layer1_3340366hid_layer1_3340368*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer1_layer_call_and_return_conditional_losses_3340185�
"hid_layer2/StatefulPartitionedCallStatefulPartitionedCall+hid_layer1/StatefulPartitionedCall:output:0hid_layer2_3340371hid_layer2_3340373*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_hid_layer2_layer_call_and_return_conditional_losses_3340202�
!out_layer/StatefulPartitionedCallStatefulPartitionedCall+hid_layer2/StatefulPartitionedCall:output:0out_layer_3340376out_layer_3340378*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *O
fJRH
F__inference_out_layer_layer_call_and_return_conditional_losses_3340219y
IdentityIdentity*out_layer/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp#^hid_layer1/StatefulPartitionedCall#^hid_layer2/StatefulPartitionedCall"^out_layer/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������: : : : : : 2H
"hid_layer1/StatefulPartitionedCall"hid_layer1/StatefulPartitionedCall2H
"hid_layer2/StatefulPartitionedCall"hid_layer2/StatefulPartitionedCall2F
!out_layer/StatefulPartitionedCall!out_layer/StatefulPartitionedCall:Y U
'
_output_shapes
:���������
*
_user_specified_namehid_layer1_input"�L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
M
hid_layer1_input9
"serving_default_hid_layer1_input:0���������=
	out_layer0
StatefulPartitionedCall:0���������tensorflow/serving/predict:�N
�
layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
	optimizer
	variables
trainable_variables
regularization_losses
	keras_api
	__call__
*
&call_and_return_all_conditional_losses
_default_save_signature

signatures"
_tf_keras_sequential
�

kernel
bias
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses"
_tf_keras_layer
�

kernel
bias
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses"
_tf_keras_layer
�

kernel
bias
	variables
 trainable_variables
!regularization_losses
"	keras_api
#__call__
*$&call_and_return_all_conditional_losses"
_tf_keras_layer
�
%iter

&beta_1

'beta_2
	(decay
)learning_ratemDmEmFmGmHmIvJvKvLvMvNvO"
	optimizer
J
0
1
2
3
4
5"
trackable_list_wrapper
J
0
1
2
3
4
5"
trackable_list_wrapper
 "
trackable_list_wrapper
�
*non_trainable_variables

+layers
,metrics
-layer_regularization_losses
.layer_metrics
	variables
trainable_variables
regularization_losses
	__call__
_default_save_signature
*
&call_and_return_all_conditional_losses
&
"call_and_return_conditional_losses"
_generic_user_object
�2�
/__inference_sequential_17_layer_call_fn_3340241
/__inference_sequential_17_layer_call_fn_3340405
/__inference_sequential_17_layer_call_fn_3340422
/__inference_sequential_17_layer_call_fn_3340342�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340448
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340474
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340362
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340382�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�B�
"__inference__wrapped_model_3340166hid_layer1_input"�
���
FullArgSpec
args� 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
,
/serving_default"
signature_map
$:"	�2hid_layer1/kernel
:�2hid_layer1/bias
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
0non_trainable_variables

1layers
2metrics
3layer_regularization_losses
4layer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
�2�
,__inference_hid_layer1_layer_call_fn_3340502�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
G__inference_hid_layer1_layer_call_and_return_conditional_losses_3340513�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
%:#
��2hid_layer2/kernel
:�2hid_layer2/bias
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
5non_trainable_variables

6layers
7metrics
8layer_regularization_losses
9layer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
�2�
,__inference_hid_layer2_layer_call_fn_3340522�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
G__inference_hid_layer2_layer_call_and_return_conditional_losses_3340533�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
#:!	�2out_layer/kernel
:2out_layer/bias
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
:non_trainable_variables

;layers
<metrics
=layer_regularization_losses
>layer_metrics
	variables
 trainable_variables
!regularization_losses
#__call__
*$&call_and_return_all_conditional_losses
&$"call_and_return_conditional_losses"
_generic_user_object
�2�
+__inference_out_layer_layer_call_fn_3340542�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
F__inference_out_layer_layer_call_and_return_conditional_losses_3340553�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
 "
trackable_list_wrapper
5
0
1
2"
trackable_list_wrapper
'
?0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
%__inference_signature_wrapper_3340493hid_layer1_input"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
N
	@total
	Acount
B	variables
C	keras_api"
_tf_keras_metric
:  (2total
:  (2count
.
@0
A1"
trackable_list_wrapper
-
B	variables"
_generic_user_object
):'	�2Adam/hid_layer1/kernel/m
#:!�2Adam/hid_layer1/bias/m
*:(
��2Adam/hid_layer2/kernel/m
#:!�2Adam/hid_layer2/bias/m
(:&	�2Adam/out_layer/kernel/m
!:2Adam/out_layer/bias/m
):'	�2Adam/hid_layer1/kernel/v
#:!�2Adam/hid_layer1/bias/v
*:(
��2Adam/hid_layer2/kernel/v
#:!�2Adam/hid_layer2/bias/v
(:&	�2Adam/out_layer/kernel/v
!:2Adam/out_layer/bias/v�
"__inference__wrapped_model_3340166z9�6
/�,
*�'
hid_layer1_input���������
� "5�2
0
	out_layer#� 
	out_layer����������
G__inference_hid_layer1_layer_call_and_return_conditional_losses_3340513]/�,
%�"
 �
inputs���������
� "&�#
�
0����������
� �
,__inference_hid_layer1_layer_call_fn_3340502P/�,
%�"
 �
inputs���������
� "������������
G__inference_hid_layer2_layer_call_and_return_conditional_losses_3340533^0�-
&�#
!�
inputs����������
� "&�#
�
0����������
� �
,__inference_hid_layer2_layer_call_fn_3340522Q0�-
&�#
!�
inputs����������
� "������������
F__inference_out_layer_layer_call_and_return_conditional_losses_3340553]0�-
&�#
!�
inputs����������
� "%�"
�
0���������
� 
+__inference_out_layer_layer_call_fn_3340542P0�-
&�#
!�
inputs����������
� "�����������
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340362rA�>
7�4
*�'
hid_layer1_input���������
p 

 
� "%�"
�
0���������
� �
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340382rA�>
7�4
*�'
hid_layer1_input���������
p

 
� "%�"
�
0���������
� �
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340448h7�4
-�*
 �
inputs���������
p 

 
� "%�"
�
0���������
� �
J__inference_sequential_17_layer_call_and_return_conditional_losses_3340474h7�4
-�*
 �
inputs���������
p

 
� "%�"
�
0���������
� �
/__inference_sequential_17_layer_call_fn_3340241eA�>
7�4
*�'
hid_layer1_input���������
p 

 
� "�����������
/__inference_sequential_17_layer_call_fn_3340342eA�>
7�4
*�'
hid_layer1_input���������
p

 
� "�����������
/__inference_sequential_17_layer_call_fn_3340405[7�4
-�*
 �
inputs���������
p 

 
� "�����������
/__inference_sequential_17_layer_call_fn_3340422[7�4
-�*
 �
inputs���������
p

 
� "�����������
%__inference_signature_wrapper_3340493�M�J
� 
C�@
>
hid_layer1_input*�'
hid_layer1_input���������"5�2
0
	out_layer#� 
	out_layer���������