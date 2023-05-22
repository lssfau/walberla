{% extends "PackingKernels.tmpl.h" %}

{% block AdditionalPublicDeclarations %}
{% if target is equalto 'cpu' -%}
   using MaskField_T = GhostLayerField< uint32_t, 1 >;
{%- elif target is equalto 'gpu' -%}
   using MaskField_T = gpu::GPUField< uint32_t >;
{%- endif %}


   /**
   * Unpacks and uniformly redistributes populations coming from a coarse block onto
   * the fine grid.
   */
   void unpackRedistribute(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
      "unsigned char * inBuffer", kernels['unpackRedistribute'].kernel_selection_parameters,
       ["gpuStream_t stream = nullptr"] if is_gpu else []]
      | type_identifier_list -}}
   ) const;

   /**
   * Partially coalesces and packs populations streaming from a fine block into a coarse block
   */
   void packPartialCoalescence(
      {{- [ "PdfField_T * " + src_field.name, "MaskField_T * " + mask_field.name, "CellInterval & ci",
      "unsigned char * outBuffer", kernels['packPartialCoalescence'].kernel_selection_parameters,
          ["gpuStream_t stream = nullptr"] if is_gpu else []]
      | type_identifier_list -}}
   ) const;

   /**
    * Prepares a coarse block for coalescence by setting every population that must be coalesced from fine blocks
    * to zero.
    */
   void zeroCoalescenceRegion(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
      kernels['zeroCoalescenceRegion'].kernel_selection_parameters,
          ["gpuStream_t stream = nullptr"] if is_gpu else []]
      | type_identifier_list -}}
   ) const;

   /**
   * Unpacks and coalesces populations coming from a fine block onto the fine grid
   */
   void unpackCoalescence(
      {{- [ "PdfField_T * " + dst_field.name, "CellInterval & ci",
      "unsigned char * inBuffer", kernels['unpackCoalescence'].kernel_selection_parameters,
          ["gpuStream_t stream = nullptr"] if is_gpu else []]
      | type_identifier_list -}}
   ) const;

   /**
    * Returns the number of bytes that will be unpacked to the cell interval
    * when using unpackRedistribute. This is 2^{-d} of the data that would be
    * unpacked during same-level communication.
    * @param ci  The cell interval
    * @return    The required size of the buffer, in bytes
    */
   uint_t redistributeSize(CellInterval & ci) const {
      return size(ci) >> {{dimension}};
   }

   /**
    * Returns the number of bytes that will be packed from the cell interval
    * when using packPartialCoalescence.
    * @param ci  The cell interval
    * @param dir The communication direction
    * @return    The required size of the buffer, in bytes
    */
   uint_t partialCoalescenceSize(CellInterval & ci, stencil::Direction dir) const {
      return size(ci, dir) >> {{dimension}};
   }
{% endblock %}
