{% extends "PackingKernels.tmpl.cpp" %}

{% block AdditionalKernelDefinitions %}
{{ kernels['unpackRedistribute']    | generate_definitions }}
{{ kernels['packPartialCoalescence']    | generate_definitions }}
{{ kernels['zeroCoalescenceRegion']    | generate_definitions }}
{{ kernels['unpackCoalescence']    | generate_definitions }}
{% endblock %}

{% block AdditionalDefinitions %}

void {{class_name}}::unpackRedistribute(
   {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
   "unsigned char * inBuffer", kernels['unpackDirection'].kernel_selection_parameters,
       ["gpuStream_t stream"] if is_gpu else []]
   | type_identifier_list -}}
) const {
   {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(inBuffer);

   {{kernels['unpackRedistribute'] | generate_call(cell_interval='ci', stream='stream') | indent(3) }}
}

void {{class_name}}::packPartialCoalescence(
   {{- [ "PdfField_T * " + src_field.name, "MaskField_T * " + mask_field.name, "CellInterval & ci",
   "unsigned char * outBuffer", kernels['packPartialCoalescence'].kernel_selection_parameters,
       ["gpuStream_t stream"] if is_gpu else []]
   | type_identifier_list -}}
) const {
   {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(outBuffer);

   {{kernels['packPartialCoalescence'] | generate_call(cell_interval='ci', stream='stream') | indent(3) }}
}

void {{class_name}}::zeroCoalescenceRegion(
   {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
   kernels['zeroCoalescenceRegion'].kernel_selection_parameters,
       ["gpuStream_t stream"] if is_gpu else []]
   | type_identifier_list -}}
) const {
   {{kernels['zeroCoalescenceRegion'] | generate_call(cell_interval='ci', stream='stream') | indent(3) }}
}

void {{class_name}}::unpackCoalescence(
   {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
   "unsigned char * inBuffer", kernels['unpackCoalescence'].kernel_selection_parameters,
       ["gpuStream_t stream"] if is_gpu else []]
   | type_identifier_list -}}
) const {
   {{dtype}} * buffer = reinterpret_cast<{{dtype}}*>(inBuffer);

   {{kernels['unpackCoalescence'] | generate_call(cell_interval='ci', stream='stream') | indent(3) }}
}

{% endblock %}
